// Author: Armin Töpfer

#include <atomic>
#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>
#include <vector>

#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>

#include <pbbam/DataSet.h>
#include <pbbam/DataSetTypes.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>

#include <Pbmm2Version.h>

#include "AlignSettings.h"
#include "MM2Helper.h"

#include "AlignWorkflow.h"

namespace PacBio {
namespace minimap2 {
namespace {
struct Summary
{
    int32_t NumAlns = 0;
    int64_t Bases = 0;
    double Similarity = 0;
};

std::tuple<std::string, std::string, std::string> CheckPositionalArgs(
    const std::vector<std::string>& args, AlignSettings* settings)
{
    if (args.size() < 2) {
        PBLOG_FATAL << "Please provide at least the input arguments: input reference output!";
        PBLOG_FATAL << "EXAMPLE: pbmm2 input.subreads.bam reference.fasta output.bam";
        std::exit(EXIT_FAILURE);
    }

    auto inputFile = args[0];
    if (!Utility::FileExists(inputFile)) {
        PBLOG_FATAL << "Input data file does not exist: " << inputFile;
        std::exit(EXIT_FAILURE);
    }
    bool isFromJson = Utility::FileExtension(inputFile) == "json";
    if (isFromJson) {
        std::ifstream ifs(inputFile);
        JSON::Json j;
        ifs >> j;
        const auto panic = [](const std::string& error) {
            PBLOG_FATAL << "JSON Datastore: " << error;
            std::exit(EXIT_FAILURE);
        };
        if (j.empty()) panic("Empty file!");
        if (j.count("files") == 0) panic("Could not find files element!");
        if (j.count("files") > 1) panic("More than ONE files element!");
        if (j["files"].empty()) panic("files element is empty!");
        if (j["files"].size() > 1) panic("files element contains more than ONE entry!");
        for (const auto& file : j["files"]) {
            if (file.count("path") == 0) panic("Could not find path element!");
            inputFile = file["path"].get<std::string>();
        }
    }
    BAM::DataSet dsInput(inputFile);
    switch (dsInput.Type()) {
        case BAM::DataSet::TypeEnum::SUBREAD: {
            if (isFromJson) {
                settings->AlignMode = AlignmentMode::SUBREADS;
                PBLOG_INFO << "Setting to SUBREAD preset";
            }
            break;
        }
        case BAM::DataSet::TypeEnum::CONSENSUS_READ: {
            if (isFromJson) {
                settings->AlignMode = AlignmentMode::CCS;
                PBLOG_INFO << "Setting to CCS preset";
            }
            break;
        }
        case BAM::DataSet::TypeEnum::TRANSCRIPT: {
            if (isFromJson) {
                settings->AlignMode = AlignmentMode::ISOSEQ;
                PBLOG_INFO << "Setting to ISOSEQ preset";
            }
            break;
        } break;
        case BAM::DataSet::TypeEnum::ALIGNMENT:
        case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
        case BAM::DataSet::TypeEnum::BARCODE:
        case BAM::DataSet::TypeEnum::REFERENCE:
        default:
            PBLOG_FATAL << "Unsupported input data file " << inputFile << " of type "
                        << BAM::DataSet::TypeToName(dsInput.Type());
            std::exit(EXIT_FAILURE);
    }

    const auto& referenceFiles = args[1];
    if (!Utility::FileExists(referenceFiles)) {
        PBLOG_FATAL << "Input reference file does not exist: " << referenceFiles;
        std::exit(EXIT_FAILURE);
    }
    std::string reference;
    if (boost::algorithm::ends_with(referenceFiles, ".mmi")) {
        reference = referenceFiles;
        PBLOG_INFO
            << "Reference input is an index file. Index parameter override options are disabled!";
    } else {
        BAM::DataSet dsRef(referenceFiles);
        switch (dsRef.Type()) {
            case BAM::DataSet::TypeEnum::REFERENCE:
                break;
            case BAM::DataSet::TypeEnum::BARCODE:
            case BAM::DataSet::TypeEnum::SUBREAD:
            case BAM::DataSet::TypeEnum::ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_READ:
            default:
                PBLOG_FATAL << "ERROR: Unsupported reference input file " << referenceFiles
                            << " of type " << BAM::DataSet::TypeToName(dsRef.Type());
                std::exit(EXIT_FAILURE);
        }
        const auto fastaFiles = dsRef.FastaFiles();
        if (fastaFiles.size() != 1) {
            PBLOG_FATAL << "Only one reference sequence allowed!";
            std::exit(EXIT_FAILURE);
        }
        reference = fastaFiles.front();
    }

    std::string out;
    if (args.size() == 3)
        out = args[2];
    else
        out = "-";

    return std::tuple<std::string, std::string, std::string>{inputFile, reference, out};
}

std::unique_ptr<BAM::internal::IQuery> BamQuery(const BAM::DataSet& ds)
{
    const auto filter = BAM::PbiFilter::FromDataSet(ds);
    std::unique_ptr<BAM::internal::IQuery> query(nullptr);
    if (filter.IsEmpty())
        query = std::make_unique<BAM::EntireFileQuery>(ds);
    else
        query = std::make_unique<BAM::PbiFilterQuery>(filter, ds);
    return query;
}

void CreateDataSet(const BAM::DataSet& originalInputDataset, const std::string& outputFile)
{
    using BAM::DataSet;
    std::string metatype;
    std::string outputType;
    const auto type = originalInputDataset.Type();
    switch (type) {
        case BAM::DataSet::TypeEnum::SUBREAD:
            metatype = "PacBio.AlignmentFile.AlignmentBamFile";
            outputType = "alignmentset";
            break;
        case BAM::DataSet::TypeEnum::CONSENSUS_READ:
            metatype = "PacBio.AlignmentFile.ConsensusAlignmentBamFile";
            outputType = "consensusalignmentset";
            break;
        default:
            throw std::runtime_error("Unsupported input type");
    }
    DataSet ds(BAM::DataSet::TypeEnum::ALIGNMENT);
    ds.Attribute("xmlns:pbdm") = "http://pacificbiosciences.com/PacBioDataModel.xsd";
    ds.Attribute("xmlns:pbmeta") = "http://pacificbiosciences.com/PacBioCollectionMetadata.xsd";
    ds.Attribute("xmlns:pbpn") = "http://pacificbiosciences.com/PacBioPartNumbers.xsd";
    ds.Attribute("xmlns:pbrk") = "http://pacificbiosciences.com/PacBioReagentKit.xsd";
    ds.Attribute("xmlns:pbsample") = "http://pacificbiosciences.com/PacBioSampleInfo.xsd";
    ds.Attribute("xmlns:pbbase") = "http://pacificbiosciences.com/PacBioBaseDataModel.xsd";

    std::string fileName = outputFile;
    if (fileName.find("/") != std::string::npos) {
        std::vector<std::string> splits;
        boost::split(splits, fileName, boost::is_any_of("/"));
        fileName = splits.back();
    }
    BAM::ExternalResource resource(metatype, fileName + ".bam");
    ds.ExternalResources().Add(resource);
    ds.Name(fileName);
    ds.TimeStampedName(fileName + "-" + BAM::CurrentTimestamp());
    std::ofstream dsOut(outputFile + "." + outputType + ".xml");
    ds.SaveToStream(dsOut);
}

std::string OutputFilePrefix(const std::string& outputFile)
{
    // Check if output type is a dataset
    const std::string outputExt = Utility::FileExtension(outputFile);
    std::string prefix = outputFile;
    if (outputExt == "xml") {
        boost::ireplace_last(prefix, ".consensusalignmentset.xml", "");
        boost::ireplace_last(prefix, ".alignmentset.xml", "");
    } else if (outputExt == "bam") {
        boost::ireplace_last(prefix, ".bam", "");
        boost::ireplace_last(prefix, ".subreads", "");
    } else {
        PBLOG_FATAL << "Unknown file extension for output file: " << outputFile;
        std::exit(EXIT_FAILURE);
    }
    return prefix;
}
}  // namespace

int AlignWorkflow::Runner(const CLI::Results& options)
{
    std::ofstream logStream;
    const Logging::LogLevel logLevel(options.IsFromRTC() ? options.LogLevel()
                                                         : options["log_level"].get<std::string>());
    {
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
    }

    AlignSettings settings(options);

    BAM::DataSet qryFile;
    std::string refFile;
    std::string outFile;
    std::string alnFile{"-"};

    std::tie(qryFile, refFile, outFile) =
        CheckPositionalArgs(options.PositionalArguments(), &settings);

    std::string outputFilePrefix;
    bool outputIsXML = false;
    if (outFile != "-") {
        outputFilePrefix = OutputFilePrefix(outFile);
        outputIsXML = Utility::FileExtension(outFile) == "xml";
        if (outputIsXML)
            alnFile = outputFilePrefix + ".bam";
        else
            alnFile = outFile;

        if (Utility::FileExists(alnFile))
            PBLOG_WARN << "Warning: Overwriting existing output file: " << alnFile;
        if (alnFile != outFile && Utility::FileExists(outFile))
            PBLOG_WARN << "Warning: Overwriting existing output file: " << outFile;
    }

    const FilterFunc filter = [&settings](const AlignedRecord& aln) {
        if (aln.Span <= 0 || aln.Span < settings.MinAlignmentLength) return false;
        if (settings.MinAccuracy <= 0) return true;
        if (aln.Similarity < settings.MinAccuracy) return false;
        return true;
    };

    MM2Helper mm2helper(refFile, settings);

    Summary s;
    int64_t alignedReads = 0;

    {
        auto qryRdr = BamQuery(qryFile);

        const auto bamFiles = qryFile.BamFiles();
        auto hdr = bamFiles.front().Header();
        for (size_t i = 1; i < bamFiles.size(); ++i)
            hdr += bamFiles.at(i).Header();

        if (!settings.SampleName.empty()) {
            auto rgs = hdr.ReadGroups();
            hdr.ClearReadGroups();
            for (auto& rg : rgs) {
                rg.Sample(settings.SampleName);
                hdr.AddReadGroup(rg);
            }
        }

        const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
        for (const auto& si : mm2helper.SequenceInfos())
            hdr.AddSequence(si);
        hdr.AddProgram(BAM::ProgramInfo("pbmm2").Name("pbmm2").Version(version).CommandLine(
            "pbmm2 " + settings.CLI));

        PacBio::Parallel::FireAndForget faf(settings.NumThreads, 3);
        BAM::BamWriter out(alnFile, hdr);

        int32_t i = 0;
        const int32_t chunkSize = settings.ChunkSize;
        auto records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);

        std::mutex outputMutex;
        int64_t alignedRecords = 0;
        std::atomic_int waiting{0};
        auto Submit = [&](const std::unique_ptr<std::vector<BAM::BamRecord>>& recs) {
            int32_t aligned = 0;
            auto output = mm2helper.Align(recs, filter, &aligned);
            if (output) {
                std::lock_guard<std::mutex> lock(outputMutex);
                alignedReads += aligned;
                for (const auto& aln : *output) {
                    s.Bases += aln.NumAlignedBases;
                    s.Similarity += aln.Similarity;
                    ++s.NumAlns;
                    out.Write(aln.Record);
                    if (++alignedRecords % 1000 == 0) {
                        PBLOG_DEBUG << "#Reads, #Alignments: " << alignedReads << ", "
                                    << alignedRecords;
                    }
                }
            }
            waiting--;
        };

        if (settings.MedianFilter) {
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
                if (logLevel == Logging::LogLevel::DEBUG) {
                    std::ostringstream ss;
                    for (size_t x = 0; x < tmp.size(); ++x) {
                        const auto& ra = tmp[x];
                        if (x == mid) ss << '[';
                        ss << ra.Length << (ra.FullLength ? 'F' : 'S');
                        if (x == mid) ss << ']';
                        ss << ' ';
                    }
                    PBLOG_DEBUG << "Median filter " << tmp.at(mid).Record.MovieName() << '/'
                                << tmp.at(mid).Record.HoleNumber() << ": " << ss.str();
                }
                return tmp.at(mid).Record;
            };

            std::vector<RecordAnnotated> ras;
            const auto Flush = [&]() {
                if (!ras.empty()) (*records)[i++] = PickMedian(std::move(ras));
            };
            for (auto& record : *qryRdr) {
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
        } else {
            while (qryRdr->GetNext((*records)[i++])) {
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
    }
    PBLOG_INFO << "Number of Aligned Reads: " << alignedReads;
    PBLOG_INFO << "Number of Alignments: " << s.NumAlns;
    PBLOG_INFO << "Number of Bases: " << s.Bases;
    PBLOG_INFO << "Mean Concordance (mapped) : "
               << std::round(1000.0 * s.Similarity / s.NumAlns) / 10.0 << "%";

    if (outputIsXML) CreateDataSet(qryFile, outputFilePrefix);

    return EXIT_SUCCESS;
}  // namespace minimap2
}  // namespace minimap2
}  // namespace PacBio
