#include "AlignWorkflow.h"

#include "AbortException.h"
#include "AlignSettings.h"
#include "BamIndex.h"
#include "InputOutputUX.h"
#include "SampleNames.h"
#include "StreamWriters.h"
#include "bam_sort.h"

#include <pbbam/BamWriter.h>
#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastqReader.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/virtual/ZmwReadStitcher.h>
#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/json/JSON.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/reports/Report.h>
#include <pbcopper/utility/FileUtils.h>
#include <pbcopper/utility/MemoryConsumption.h>
#include <pbcopper/utility/Stopwatch.h>
#include <pbmm2/MM2Helper.h>

#include <sys/stat.h>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>
#include <tuple>
#include <vector>

namespace PacBio {
namespace minimap2 {

int AlignWorkflow::Runner(const CLI_v2::Results& options)
{
    Utility::Stopwatch startTime;
    AlignSettings settings(options);

    UserIO uio = InputOutputUX::CheckPositionalArgs(options.PositionalArguments(), settings);

    if (!uio.isAlignedInput && (!uio.isFromXML || !uio.isFromSubreadset)) {
        if (settings.ZMW) {
            throw AbortException(
                "Option --zmw can only be used with a subreadset.xml containing subread + "
                "scraps BAM files.");
        }
        if (settings.HQRegion) {
            throw AbortException(
                "Option --hqregion can only be used with a subreadset.xml containing subread + "
                "scraps BAM files.");
        }
    }

    if (uio.isFromMmi && settings.CompressSequenceHomopolymers) {
        throw AbortException("Cannot combine --collapse-homopolymers with MMI input.");
    }

    if (uio.isFromFofn && settings.SplitBySample) {
        throw AbortException("Cannot combine --split-by-sample with fofn input.");
    }

    if (!uio.isFastaInput && !uio.isFastqInput && !settings.Rg.empty()) {
        throw AbortException("Cannot override read groups with BAM input. Remove option --rg.");
    }

    const FilterFunc filter = [&settings](const AlignedRecord& aln) {
        if (aln.Span <= 0 || aln.Span < settings.MinAlignmentLength) return false;
        if (settings.MinPercIdentity <= 0 && settings.MinPercIdentityGapComp <= 0 &&
            settings.MinPercConcordance <= 0)
            return true;
        return (aln.Identity >= settings.MinPercIdentity &&
                aln.IdentityGapComp >= settings.MinPercIdentityGapComp &&
                aln.Concordance >= settings.MinPercConcordance);
    };

    const auto CompressHomopolymers = [&](const std::string& bases) -> std::string {
        const int32_t numBases = bases.size();
        std::string compressed{bases.at(0)};
        for (int32_t i = 1; i < numBases; ++i)
            if (std::toupper(bases.at(i)) != std::toupper(bases.at(i - 1)))
                compressed += bases.at(i);
        return compressed;
    };

    Utility::Stopwatch indexTime;
    std::unique_ptr<MM2Helper> mm2helper;
    if (settings.CompressSequenceHomopolymers) {
        std::vector<BAM::FastaSequence> refs = BAM::FastaReader::ReadAll(uio.refFile);
        std::ofstream compressedRef(uio.outPrefix + ".ref.collapsed.fasta");
        for (size_t i = 0; i < refs.size(); ++i) {
            const std::string rle = CompressHomopolymers(refs[i].Bases());
            compressedRef << '>' << refs[i].Name() << '\n' << rle << '\n';
            refs[i] = BAM::FastaSequence(refs[i].Name(), std::move(rle));
        }
        mm2helper = std::make_unique<MM2Helper>(refs, settings);
    } else {
        mm2helper = std::make_unique<MM2Helper>(uio.refFile, settings);
    }
    indexTime.Freeze();
    Utility::Stopwatch alignmentTime;

    Summary s;
    int64_t alignedReads = 0;

    const auto BamQueryFile = [](const std::string& file) {
        try {
            const auto filter = BAM::PbiFilter::FromDataSet(file);
            std::unique_ptr<BAM::internal::IQuery> query(nullptr);
            if (filter.IsEmpty())
                query = std::make_unique<BAM::EntireFileQuery>(file);
            else
                query = std::make_unique<BAM::PbiFilterQuery>(filter, file);
            return query;
        } catch (...) {
            throw AbortException(UNKNOWN_FILE_TYPES);
        }
    };

    std::unique_ptr<StreamWriters> writers;
    {
        static const std::string fallbackSampleName{"UnnamedSample"};

        MovieToSampleToInfix mtsti;
        if (!uio.isFastaInput && !uio.isFastqInput)
            mtsti = SampleNames::DetermineMovieToSampleToInfix(uio);

        std::map<std::string, std::set<std::string>> infixToSamples;
        for (const auto& movie_sampleInfix : mtsti)
            infixToSamples[movie_sampleInfix.second.second].insert(movie_sampleInfix.second.first);

        for (const auto& infix_samples : infixToSamples) {
            auto& infix = infix_samples.first;
            if (infix_samples.second.size() > 1) {
                int infixCounter = 0;
                for (const auto& sample : infix_samples.second) {
                    std::string newInfix = infix + '-' + std::to_string(infixCounter);
                    for (auto& movie_sampleInfix : mtsti) {
                        if (movie_sampleInfix.second.first == sample) {
                            movie_sampleInfix.second.second = newInfix;
                        }
                    }
                    ++infixCounter;
                }
            }
        }

        std::string fastxRgId = "default";
        BAM::BamHeader hdr = SampleNames::GenerateBamHeader(settings, uio, mtsti, fastxRgId);
        for (const auto& si : mm2helper->SequenceInfos())
            hdr.AddSequence(si);

        PacBio::Parallel::FireAndForget faf(settings.NumThreads, 3);

        writers = std::make_unique<StreamWriters>(
            hdr, uio.outPrefix, settings.SplitBySample, settings.Sort, settings.BamIdx,
            settings.SortThreads, settings.NumThreads, settings.SortMemory);

        int32_t i = 0;
        const int32_t chunkSize = settings.ChunkSize;
        auto records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);

        std::mutex outputMutex;
        int64_t alignedRecords = 0;
        std::atomic_int waiting{0};
        const auto firstTime = std::chrono::steady_clock::now();
        auto lastTime = std::chrono::steady_clock::now();
        auto Submit = [&](const std::unique_ptr<std::vector<BAM::BamRecord>>& recs) {
            const auto Strip = [](BAM::BamRecord& record) {
                constexpr std::array<const char[3], 25> STRIP_TAGS = {
                    "dq", "dt", "fi", "fn", "fp", "ip", "iq", "mq", "pa", "pc", "pd", "pe", "pg",
                    "pm", "pq", "pt", "pv", "pw", "px", "ri", "rn", "rp", "sf", "sq", "st"};
                auto& impl = record.Impl();
                for (const auto& tag : STRIP_TAGS) {
                    impl.RemoveTag(tag);
                }
            };
            if (settings.Strip) {
                for (auto& r : *recs)
                    Strip(r);
            }
            if (settings.CompressSequenceHomopolymers) {
                const auto Compress = [&](BAM::BamRecord& record) {
                    std::string newSeq = CompressHomopolymers(record.Sequence());
                    record.Impl().SetSequenceAndQualities(newSeq);
                    if (record.HasQueryStart() && record.HasQueryEnd()) {
                        record.Impl().EditTag(
                            "qe", record.QueryStart() + static_cast<int32_t>(newSeq.size()));
                    }
                };
                for (auto& r : *recs) {
                    Compress(r);
                    Strip(r);
                }
            }
            int32_t aligned = 0;
            try {
                auto output = mm2helper->Align(recs, filter, &aligned);
                if (output) {
                    std::lock_guard<std::mutex> lock(outputMutex);
                    alignedReads += aligned;
                    for (auto& aln : *output) {
                        if (!settings.OutputUnmapped && !aln.IsAligned) continue;
                        if (aln.IsAligned) {
                            s.Lengths.emplace_back(aln.NumAlignedBases);
                            s.Bases += aln.NumAlignedBases;
                            s.Concordance += aln.Concordance;
                            s.Identity += aln.Identity;
                            s.IdentityGapComp += aln.IdentityGapComp;
                            ++s.NumAlns;
                            if (settings.MinPercConcordance <= 0) aln.Record.Impl().RemoveTag("mc");
                            if (settings.MinPercIdentityGapComp <= 0)
                                aln.Record.Impl().RemoveTag("mg");
                            if (settings.MinPercIdentity <= 0) aln.Record.Impl().RemoveTag("mi");
                        } else {
                            aln.Record.Impl().MapQuality(0);
                        }
                        const std::string movieName = aln.Record.MovieName();
                        const auto& sampleInfix = mtsti[movieName];
                        writers->at(sampleInfix.second, sampleInfix.first).Write(aln.Record);
                        if (aln.IsAligned && ++alignedRecords % settings.ChunkSize == 0) {
                            const auto now = std::chrono::steady_clock::now();
                            auto elapsedSecs =
                                std::chrono::duration_cast<std::chrono::seconds>(now - lastTime)
                                    .count();
                            if (elapsedSecs > 5) {
                                lastTime = now;
                                auto elapsedSecTotal =
                                    std::chrono::duration_cast<std::chrono::seconds>(now -
                                                                                     firstTime)
                                        .count() /
                                    60.0;
                                auto alnsPerMin = std::round(alignedReads / elapsedSecTotal);
                                PBLOG_DEBUG << "#Reads, #Aln, #RPM: " << alignedReads << ", "
                                            << alignedRecords << ", " << alnsPerMin;
                            }
                        }
                    }
                }
            } catch (...) {
                std::cerr << "ERROR" << std::endl;
            }
            waiting--;
        };

        const auto FastxToUnalignedBam = [&hdr, &fastxRgId](const std::string& seq,
                                                            const std::string& name,
                                                            const std::string& qual) {
            BAM::BamRecord record(hdr);
            record.Impl().SetSequenceAndQualities(seq, qual);
            record.ReadGroupId(fastxRgId);
            record.Impl().Name(name);
            record.Impl().AddTag("qs", 0);
            record.Impl().AddTag("qe", static_cast<int32_t>(seq.size()));
            return record;
        };

        if (uio.isFastaInput) {
            for (const auto& f : uio.inputFiles) {
                BAM::FastaReader reader(f);
                BAM::FastaSequence fa;
                while (reader.GetNext(fa)) {
                    (*records)[i++] = FastxToUnalignedBam(fa.Bases(), fa.Name(), "");
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                }
            }
        } else if (uio.isFastqInput) {
            for (const auto& f : uio.inputFiles) {
                BAM::FastqReader reader(f);
                BAM::FastqSequence fq;
                while (reader.GetNext(fq)) {
                    (*records)[i++] =
                        FastxToUnalignedBam(fq.Bases(), fq.Name(), fq.Qualities().Fastq());
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                }
            }
        } else if (uio.isAlignedInput) {
            if (settings.MedianFilter)
                PBLOG_WARN << "Option --median-filter is ignored with aligned input!";
            if (settings.ZMW) PBLOG_WARN << "Option --zmw is ignored with aligned input!";
            if (settings.HQRegion) PBLOG_WARN << "Option --hqregion is ignored with aligned input!";

            const auto Fill = [&](const std::string& f) {
                auto reader = BamQueryFile(f);
                BAM::BamRecord tmp;
                while (reader->GetNext(tmp)) {
                    if (tmp.Impl().IsSupplementaryAlignment()) continue;
                    (*records)[i++] = std::move(tmp);
                    tmp = BAM::BamRecord();
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                }
            };
            if (uio.isFromJson) {
                Fill(uio.unpackedFromJson);
            } else {
                for (const auto& f : uio.inputFiles)
                    Fill(f);
            }
        } else if (settings.MedianFilter) {
            struct RecordAnnotated
            {
                RecordAnnotated(BAM::BamRecord record)
                    : Record(std::move(record)), Length(Record.Sequence().size())
                {
                    if (Record.HasLocalContextFlags()) {
                        const auto flags = Record.LocalContextFlags();
                        if (flags & Data::ADAPTER_BEFORE && flags & Data::ADAPTER_AFTER) {
                            FullLength = true;
                        }
                    }
                }

                BAM::BamRecord Record;
                int32_t Length;
                bool FullLength = false;
            };

            std::string movieName;
            int32_t holeNumber = -1;

            const auto PickMedian = [&](std::vector<RecordAnnotated> tmp) {
                bool hasFullLength = false;
                for (const auto& ra : tmp) {
                    if (ra.FullLength) {
                        hasFullLength = true;
                        break;
                    }
                }
                if (hasFullLength) {
                    std::vector<RecordAnnotated> fullLengths;
                    fullLengths.reserve(tmp.size());
                    for (auto&& ra : tmp)
                        if (ra.FullLength) fullLengths.emplace_back(std::move(ra));
                    tmp.swap(fullLengths);
                }

                std::stable_sort(tmp.begin(), tmp.end(),
                                 [&](const RecordAnnotated& l, const RecordAnnotated& r) {
                                     return std::tie(l.FullLength, l.Length) <
                                            std::tie(r.FullLength, r.Length);
                                 });
                size_t mid = tmp.size() / 2;
                if (options.LogLevel() == Logging::LogLevel::TRACE) {
                    std::ostringstream ss;
                    for (size_t x = 0; x < tmp.size(); ++x) {
                        const auto& ra = tmp[x];
                        if (x == mid) ss << '[';
                        ss << ra.Length << (ra.FullLength ? 'F' : 'S');
                        if (x == mid) ss << ']';
                        ss << ' ';
                    }
                    PBLOG_TRACE << "Median filter " << tmp.at(mid).Record.MovieName() << '/'
                                << tmp.at(mid).Record.HoleNumber() << ": " << ss.str();
                }
                return tmp.at(mid).Record;
            };

            std::vector<RecordAnnotated> ras;
            const auto Flush = [&]() {
                if (!ras.empty()) (*records)[i++] = PickMedian(std::move(ras));
            };

            const auto Fill = [&](const std::string& f) {
                auto reader = BamQueryFile(f);
                for (auto& record : *reader) {
                    const auto nextHoleNumber = record.HoleNumber();
                    const auto nextMovieName = record.MovieName();
                    if (holeNumber != nextHoleNumber || movieName != nextMovieName) {
                        Flush();
                        holeNumber = nextHoleNumber;
                        movieName = nextMovieName;
                        ras = std::vector<RecordAnnotated>();
                    }
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                    ras.emplace_back(record);
                }
                Flush();
            };
            if (uio.isFromJson) {
                Fill(uio.unpackedFromJson);
            } else {
                for (const auto& f : uio.inputFiles)
                    Fill(f);
            }
        } else if (settings.HQRegion) {
            const auto Fill = [&](const std::string& f) {
                BAM::ZmwReadStitcher reader(f);
                while (reader.HasNext()) {
                    auto r = reader.Next();
                    if (r.HasVirtualRegionType(BAM::VirtualRegionType::HQREGION)) {
                        auto hqs = r.VirtualRegionsTable(BAM::VirtualRegionType::HQREGION);
                        if (hqs.empty()) {
                            PBLOG_WARN << "Skipping ZMW record " << r.FullName()
                                       << " missing HQ region";
                        } else if (hqs.size() > 1) {
                            PBLOG_WARN << "ZMW record " << r.FullName()
                                       << " has more than one HQ region, will use first";
                        } else if (hqs.at(0).beginPos == hqs.at(0).endPos) {
                            PBLOG_WARN << "ZMW record " << r.FullName()
                                       << " has a zero-length HQ region";
                            continue;
                        }
                        PBLOG_DEBUG << "clipping " << r.FullName() << " from " << hqs.at(0).beginPos
                                    << " to " << hqs.at(0).endPos;
                        r.Clip(BAM::ClipType::CLIP_TO_QUERY, hqs.at(0).beginPos, hqs.at(0).endPos);
                    }
                    (*records)[i++] = std::move(r);
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                }
            };
            if (uio.isFromJson) {
                Fill(uio.unpackedFromJson);
            } else {
                for (const auto& f : uio.inputFiles)
                    Fill(f);
            }
        } else if (settings.ZMW) {
            const auto Fill = [&](const std::string& f) {
                BAM::ZmwReadStitcher reader(f);
                while (reader.HasNext()) {
                    (*records)[i++] = reader.Next();
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                }
            };
            if (uio.isFromJson) {
                Fill(uio.unpackedFromJson);
            } else {
                for (const auto& f : uio.inputFiles)
                    Fill(f);
            }
        } else {
            const auto Fill = [&](const std::string& f) {
                auto reader = BamQueryFile(f);
                while (reader->GetNext((*records)[i++])) {
                    if (i >= chunkSize) {
                        waiting++;
                        faf.ProduceWith(Submit, std::move(records));
                        records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                        i = 0;
                    }
                }
                if (i > 0) i--;
            };
            if (uio.isFromJson) {
                Fill(uio.unpackedFromJson);
            } else {
                for (const auto& f : uio.inputFiles)
                    Fill(f);
            }
        }
        // terminal records, if they exist
        if (i > 0) {
            records->resize(i);
            waiting++;
            faf.ProduceWith(Submit, std::move(records));
        }

        faf.Finalize();

        while (waiting) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        if (settings.Sort)
            PBLOG_DEBUG << "Alignment finished, merging sorted chunks using "
                        << (settings.NumThreads + settings.SortThreads) << " threads.";
    }

    alignmentTime.Freeze();
    const auto sort_baiTimings = writers->Close();

    int32_t maxMappedLength = 0;
    for (const auto& l : s.Lengths) {
        maxMappedLength = std::max(maxMappedLength, l);
    }

    const auto DenomNumAlns = std::max(1, s.NumAlns);  // Avoid dividing by zero later
    double meanMappedConcordance = 1.0 * s.Concordance / DenomNumAlns;
    double meanIdentity = 1.0 * s.Identity / DenomNumAlns;
    double meanIdentityGapComp = 1.0 * s.IdentityGapComp / DenomNumAlns;
    double meanMappedReadLength = (1.0 * s.Bases / DenomNumAlns);

    std::string pbiTiming;
    if (uio.isToXML || uio.isToJson)
        pbiTiming = writers->WriteDatasetsJson(uio, s, settings.SplitBySample);
    else if (settings.CreatePbi)
        pbiTiming = writers->ForcePbiOutput();

    PBLOG_INFO << "Mapped Reads: " << alignedReads;
    PBLOG_INFO << "Alignments: " << s.NumAlns;
    PBLOG_INFO << "Mapped Bases: " << s.Bases;
    if (settings.MinPercConcordance > 0)
        PBLOG_INFO << "Mean Mapped Concordance: " << meanMappedConcordance << "%";
    if (settings.MinPercIdentity > 0)
        PBLOG_INFO << "Mean Sequence Identity: " << meanIdentity << "%";
    if (settings.MinPercIdentityGapComp > 0)
        PBLOG_INFO << "Mean Gap-Compressed Sequence Identity: " << meanIdentityGapComp << "%";
    PBLOG_INFO << "Max Mapped Read Length: " << maxMappedLength;
    PBLOG_INFO << "Mean Mapped Read Length: " << meanMappedReadLength;

    PBLOG_INFO << "Index Build/Read Time: " << indexTime.ElapsedTime();
    PBLOG_INFO << "Alignment Time: " << alignmentTime.ElapsedTime();
    if (!sort_baiTimings.first.empty()) PBLOG_INFO << "Sort Merge Time: " << sort_baiTimings.first;
    if (!sort_baiTimings.second.empty() && settings.BamIdx != BamIndex::_from_string("NONE")) {
        PBLOG_INFO << settings.BamIdx._to_string()
                   << " Generation Time: " << sort_baiTimings.second;
    }
    if (!pbiTiming.empty()) PBLOG_INFO << "PBI Generation Time: " << pbiTiming;
    PBLOG_INFO << "Run Time: " << startTime.ElapsedTime();
    PBLOG_INFO << "CPU Time: "
               << Utility::Stopwatch::PrettyPrintNanoseconds(
                      static_cast<int64_t>(Utility::Stopwatch::CpuTime() * 1000 * 1000 * 1000));
    int64_t const peakRss = PacBio::Utility::MemoryConsumption::PeakRss();
    double const peakRssGb = peakRss / 1024.0 / 1024.0 / 1024.0;
    PBLOG_INFO << "Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb << " GB";

    if (!settings.ReportFileJson.empty()) {
        Reports::Report jsonReport{"mapping_stats", "Mapping Statistics"};
        jsonReport.AddAttribute(
            {"mapped_reads_n", Reports::ReportValue{alignedReads}, "Mapped Reads"});
        jsonReport.AddAttribute(
            {"mapped_alignments_n", Reports::ReportValue{s.NumAlns}, "Alignments"});
        jsonReport.AddAttribute({"mapped_bases_n", Reports::ReportValue{s.Bases}, "Mapped Bases"});
        jsonReport.AddAttribute(
            {"blast_identity", Reports::ReportValue{meanIdentity}, "Mean Sequence Identity"});
        jsonReport.AddAttribute({"pb_identity", Reports::ReportValue{meanMappedConcordance},
                                 "Mean Mapped Concordance"});
        jsonReport.AddAttribute({"gap_compressed_identity",
                                 Reports::ReportValue{meanIdentityGapComp},
                                 "Mean Gap-Compressed Sequence Identity"});
        jsonReport.AddAttribute({"mapped_readlength_mean",
                                 Reports::ReportValue{meanMappedReadLength},
                                 "Mean Mapped Read Length"});
        jsonReport.AddAttribute({"mapped_readlength_max", Reports::ReportValue{maxMappedLength},
                                 "Max Mapped Read Length"});
        PBLOG_INFO << "Writing mapping stats report JSON to " << settings.ReportFileJson;
        jsonReport.Print(settings.ReportFileJson);
    }

    return EXIT_SUCCESS;
}
}  // namespace minimap2
}  // namespace PacBio
