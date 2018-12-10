// Author: Armin Töpfer

#include <sys/stat.h>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <tuple>
#include <vector>

#include <pbcopper/PbcopperMakeUnique.h>
#include <pbcopper/json/JSON.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamWriter.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastqReader.h>

#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/virtual/ZmwReadStitcher.h>

#include <pbcopper/parallel/FireAndForget.h>

#include <mmpriv.h>

#include <pbmm2/MM2Helper.h>
#include "AlignSettings.h"
#include "InputOutputUX.h"
#include "SampleNames.h"
#include "StreamWriters.h"
#include "Timer.h"
#include "bam_sort.h"

#include "AlignWorkflow.h"

namespace PacBio {
namespace minimap2 {
namespace {
Logging::LogLevel CreateLogger(const CLI::Results& options, std::ofstream& logStream)
{
    const Logging::LogLevel logLevel(options.IsFromRTC() ? options.LogLevel()
                                                         : options["log_level"].get<std::string>());
    const std::string logFile = options["log_file"];

    using Logger = PacBio::Logging::Logger;

    Logger* logger;
    if (!logFile.empty()) {
        logStream.open(logFile);
        logger = &Logger::Default(new Logger(logStream, logLevel));
    } else {
        logger = &Logger::Default(new Logger(std::cerr, logLevel));
    }
    PacBio::Logging::InstallSignalHandlers(*logger);
    return logLevel;
}
}  // namespace

int AlignWorkflow::Runner(const CLI::Results& options)
{
    const Timer startTime;
    std::ofstream logStream;
    const Logging::LogLevel logLevel = CreateLogger(options, logStream);

    AlignSettings settings(options);

    UserIO uio = InputOutputUX::CheckPositionalArgs(options.PositionalArguments(), settings);
    BAM::DataSet inFile;
    try {
        if (!uio.isFastaInput && !uio.isFastqInput) inFile = BAM::DataSet(uio.inFile);
    } catch (...) {
        PBLOG_FATAL << UNKNOWN_FILE_TYPES;
        std::exit(EXIT_FAILURE);
    }

    if (!uio.isAlignedInput && (!uio.isFromXML || !uio.isFromSubreadset)) {
        if (settings.ZMW) {
            PBLOG_FATAL
                << "Option --zmw can only be used with a subreadset.xml containing subread + "
                   "scraps BAM files.";
            std::exit(EXIT_FAILURE);
        }
        if (settings.HQRegion) {
            PBLOG_FATAL
                << "Option --hqregion can only be used with a subreadset.xml containing subread + "
                   "scraps BAM files.";
            std::exit(EXIT_FAILURE);
        }
    }

    if (!uio.isFastaInput && !uio.isFastqInput && !settings.Rg.empty()) {
        PBLOG_FATAL << "Cannot override read groups with BAM input. Remove option --rg.";
        std::exit(EXIT_FAILURE);
    }

    const FilterFunc filter = [&settings](const AlignedRecord& aln) {
        if (aln.Span <= 0 || aln.Span < settings.MinAlignmentLength) return false;
        if (settings.MinPercConcordance <= 0) return true;
        if (aln.Concordance < settings.MinPercConcordance) return false;
        return true;
    };

    Timer indexTime;
    MM2Helper mm2helper(uio.refFile, settings);
    indexTime.Freeze();
    Timer alignmentTime;

    Summary s;
    int64_t alignedReads = 0;

    const auto BamQuery = [&inFile]() {
        try {
            const auto filter = BAM::PbiFilter::FromDataSet(inFile);
            std::unique_ptr<BAM::internal::IQuery> query(nullptr);
            if (filter.IsEmpty())
                query = std::make_unique<BAM::EntireFileQuery>(inFile);
            else
                query = std::make_unique<BAM::PbiFilterQuery>(filter, inFile);
            return query;
        } catch (...) {
            PBLOG_FATAL << UNKNOWN_FILE_TYPES;
            std::exit(EXIT_FAILURE);
        }
    };

    std::unique_ptr<StreamWriters> writers;
    {
        static const std::string fallbackSampleName{"UnnamedSample"};

        auto mtsti = SampleNames::DetermineMovieToSampleToInfix(inFile);

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
        BAM::BamHeader hdr =
            SampleNames::GenerateBamHeader(inFile, settings, uio, mtsti, fastxRgId);
        for (const auto& si : mm2helper.SequenceInfos())
            hdr.AddSequence(si);

        PacBio::Parallel::FireAndForget faf(settings.NumThreads, 3);

        writers = std::make_unique<StreamWriters>(
            hdr, uio.outPrefix, settings.SplitBySample, settings.Sort, !settings.NoBAI,
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
            if (settings.Strip) {
                const auto Strip = [](BAM::BamRecord& record) {
                    auto& impl = record.Impl();
                    for (const auto& t :
                         {"dq", "dt", "ip", "iq", "mq", "pa", "pc", "pd", "pe", "pg", "pm", "pq",
                          "pt", "pv", "pw", "px", "sf", "sq", "st"})
                        impl.RemoveTag(t);
                };
                for (auto& r : *recs)
                    Strip(r);
            }
            int32_t aligned = 0;
            try {
                auto output = mm2helper.Align(recs, filter, &aligned);
                if (output) {
                    std::lock_guard<std::mutex> lock(outputMutex);
                    alignedReads += aligned;
                    for (const auto& aln : *output) {
                        s.Lengths.emplace_back(aln.NumAlignedBases);
                        s.Bases += aln.NumAlignedBases;
                        s.Concordance += aln.Concordance;
                        ++s.NumAlns;
                        const std::string movieName = aln.Record.MovieName();
                        const auto& sampleInfix = mtsti[movieName];
                        writers->at(sampleInfix.second, sampleInfix.first).Write(aln.Record);
                        if (++alignedRecords % settings.ChunkSize == 0) {
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
            BAM::FastaReader reader(uio.inFile);
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
        } else if (uio.isFastqInput) {
            BAM::FastqReader reader(uio.inFile);
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
        } else if (uio.isAlignedInput) {
            if (settings.MedianFilter)
                PBLOG_WARN << "Option --median-filter is ignored with aligned input!";
            if (settings.ZMW) PBLOG_WARN << "Option --zmw is ignored with aligned input!";
            if (settings.HQRegion) PBLOG_WARN << "Option --hqregion is ignored with aligned input!";

            auto reader = BamQuery();
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
        } else if (settings.MedianFilter) {
            struct RecordAnnotated
            {
                RecordAnnotated(BAM::BamRecord record)
                    : Record(std::move(record)), Length(Record.Sequence().size())
                {
                    if (Record.HasLocalContextFlags()) {
                        const auto flags = Record.LocalContextFlags();
                        if (flags & BAM::ADAPTER_BEFORE && flags & BAM::ADAPTER_AFTER) {
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
                                     return std::tie(l.FullLength, l.Length) <=
                                            std::tie(r.FullLength, r.Length);
                                 });
                size_t mid = tmp.size() / 2;
                if (logLevel == Logging::LogLevel::TRACE) {
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
            auto reader = BamQuery();
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
        } else if (settings.HQRegion) {
            BAM::ZmwReadStitcher reader(inFile);
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
                    }
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
        } else if (settings.ZMW) {
            BAM::ZmwReadStitcher reader(inFile);
            while (reader.HasNext()) {
                (*records)[i++] = reader.Next();
                if (i >= chunkSize) {
                    waiting++;
                    faf.ProduceWith(Submit, std::move(records));
                    records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                    i = 0;
                }
            }
        } else {
            auto reader = BamQuery();
            while (reader->GetNext((*records)[i++])) {
                if (i >= chunkSize) {
                    waiting++;
                    faf.ProduceWith(Submit, std::move(records));
                    records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                    i = 0;
                }
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
    double meanMappedConcordance = 1.0 * s.Concordance / s.NumAlns;
    if (settings.IsFromRTC) {
        JSON::Json root;
        root["mapped_concordance_percentage"] = meanMappedConcordance;
        root["num_aligned_reads"] = alignedReads;
        root["num_alignments"] = s.NumAlns;
        root["num_aligned_bases"] = s.Bases;
        std::ofstream out(".pbmm2_stats.json");
        out << root.dump(2);
    }

    std::string pbiTiming;
    if (uio.isToXML || uio.isToJson)
        pbiTiming =
            writers->WriteDatasetsJson(inFile, uio.outFile, uio.refFile, uio.isFromXML,
                                       uio.isToJson, s, uio.outPrefix, settings.SplitBySample);
    else if (settings.CreatePbi)
        pbiTiming = writers->ForcePbiOutput();

    PBLOG_INFO << "Mapped Reads: " << alignedReads;
    PBLOG_INFO << "Alignments: " << s.NumAlns;
    PBLOG_INFO << "Mapped Bases: " << s.Bases;
    PBLOG_INFO << "Mean Mapped Concordance: " << meanMappedConcordance << "%";
    PBLOG_INFO << "Max Mapped Read Length: " << maxMappedLength;
    PBLOG_INFO << "Mean Mapped Read Length: " << (1.0 * s.Bases / s.NumAlns);

    PBLOG_INFO << "Index Build/Read Time: " << indexTime.ElapsedTime();
    PBLOG_INFO << "Alignment Time: " << alignmentTime.ElapsedTime();
    if (!sort_baiTimings.first.empty()) PBLOG_INFO << "Sort Merge Time: " << sort_baiTimings.first;
    if (!sort_baiTimings.second.empty() && !settings.NoBAI)
        PBLOG_INFO << "BAI Generation Time: " << sort_baiTimings.second;
    if (!pbiTiming.empty()) PBLOG_INFO << "PBI Generation Time: " << pbiTiming;
    PBLOG_INFO << "Run Time: " << startTime.ElapsedTime();
    PBLOG_INFO << "CPU Time: "
               << Timer::ElapsedTimeFromSeconds(
                      static_cast<int64_t>(cputime() * 1000 * 1000 * 1000));
    PBLOG_INFO << "Peak RSS: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

    return EXIT_SUCCESS;
}
}  // namespace minimap2
}  // namespace PacBio
