// Author: Armin Töpfer

#include <sys/stat.h>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>
#include <vector>

#include <pbcopper/json/JSON.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastqReader.h>

#include <pbbam/BamFile.h>
#include <pbbam/DataSet.h>
#include <pbbam/DataSetTypes.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFile.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/virtual/ZmwReadStitcher.h>

#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pbcopper/utility/Stopwatch.h>

#include <boost/algorithm/string.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <Pbmm2Version.h>

#include <mmpriv.h>

#include <pbmm2/MM2Helper.h>
#include "AlignSettings.h"
#include "Timer.h"
#include "bam_sort.h"

#include "AlignWorkflow.h"

namespace PacBio {
namespace minimap2 {
namespace {
static const std::string UNKNOWN_FILE_TYPES =
    "Could not determine read input type(s). Please do not mix data types, such as BAM+FASTQ. File "
    "of files may only contain BAMs or datasets.";
struct Summary
{
    int32_t NumAlns = 0;
    int64_t Bases = 0;
    double Concordance = 0;
    std::vector<int32_t> Lengths;
};

static std::string UnpackJson(const std::string& jsonInputFile)
{
    std::ifstream ifs(jsonInputFile);
    JSON::Json j;
    ifs >> j;
    std::string inputFile;
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
    return inputFile;
}

enum class InputType : int
{
    READ = 0,
    FASTX,
    MMI,
    FASTA,
    FASTQ,
    UNKNOWN
};

static InputType DetermineInputTypeFastx(const std::string& in)
{
    char firstChar = ' ';
    std::ifstream f(in);
    f.read(&firstChar, 1);
    if (firstChar == '>')
        return InputType::FASTA;
    else if (firstChar == '@')
        return InputType::FASTQ;
    else {
        PBLOG_FATAL << "Unkown file type for file " << in << " starting with character "
                    << firstChar;
        std::exit(EXIT_FAILURE);
    }
}

static InputType DetermineInputTypeApprox(std::string inputFile)
{
    if (!Utility::FileExists(inputFile)) {
        PBLOG_FATAL << "Input data file does not exist: " << inputFile;
        std::exit(EXIT_FAILURE);
    }
    if (Utility::FileExtension(inputFile) == "json") inputFile = UnpackJson(inputFile);
    if (boost::algorithm::ends_with(inputFile, ".mmi")) return InputType::MMI;
    if (boost::algorithm::ends_with(inputFile, ".fq")) return InputType::FASTX;
    if (boost::algorithm::ends_with(inputFile, ".fastq")) return InputType::FASTX;
    if (boost::algorithm::ends_with(inputFile, ".fa")) return InputType::FASTX;
    if (boost::algorithm::ends_with(inputFile, ".fasta")) return InputType::FASTX;
    BAM::DataSet dsInput;
    try {
        dsInput = BAM::DataSet(inputFile);
    } catch (...) {
        PBLOG_FATAL << UNKNOWN_FILE_TYPES;
        std::exit(EXIT_FAILURE);
    }
    switch (dsInput.Type()) {
        case BAM::DataSet::TypeEnum::ALIGNMENT:
        case BAM::DataSet::TypeEnum::SUBREAD:
        case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
        case BAM::DataSet::TypeEnum::CONSENSUS_READ:
        case BAM::DataSet::TypeEnum::TRANSCRIPT_ALIGNMENT:
        case BAM::DataSet::TypeEnum::TRANSCRIPT:
            return InputType::READ;
        case BAM::DataSet::TypeEnum::BARCODE:
        case BAM::DataSet::TypeEnum::REFERENCE:
            return InputType::FASTX;
        default:
            PBLOG_FATAL << "Unsupported input data file " << inputFile << " of type "
                        << BAM::DataSet::TypeToName(dsInput.Type());
            std::exit(EXIT_FAILURE);
    }
}

std::tuple<std::string, std::string, std::string> CheckPositionalArgs(
    const std::vector<std::string>& args, AlignSettings* settings, bool* isFromJson,
    bool* isFromXML, BAM::DataSet::TypeEnum* inputType, bool* isAlignedInput, bool* isFastaInput,
    bool* isFastqInput)
{
    if (args.size() < 2) {
        PBLOG_FATAL << "Please provide at least the input arguments: input reference output!";
        PBLOG_FATAL << "EXAMPLE: pbmm2 reference.fasta input.subreads.bam output.bam";
        std::exit(EXIT_FAILURE);
    }

    std::string inputFile;
    std::string referenceFile;
    auto file0Type = DetermineInputTypeApprox(args[0]);
    auto file1Type = DetermineInputTypeApprox(args[1]);
    if (file0Type == InputType::READ &&
        (file1Type == InputType::FASTX || file1Type == InputType::MMI)) {
        inputFile = args[0];
        referenceFile = args[1];
    } else if ((file0Type == InputType::FASTX || file0Type == InputType::MMI) &&
               file1Type == InputType::READ) {
        referenceFile = args[0];
        inputFile = args[1];
    } else if ((file0Type == InputType::FASTX || file0Type == InputType::MMI) &&
               (file1Type == InputType::FASTX || file1Type == InputType::MMI)) {
        int32_t numMmi = (file0Type == InputType::MMI) + (file1Type == InputType::MMI);
        if (numMmi == 2) {
            PBLOG_FATAL << "Both input files are of type MMI. Please check your inputs.";
            std::exit(EXIT_FAILURE);
        } else if (numMmi == 1) {
            if (file0Type == InputType::MMI) {
                referenceFile = args[0];
                inputFile = args[1];
            } else {
                referenceFile = args[1];
                inputFile = args[0];
            }
            InputType inputFileType = DetermineInputTypeFastx(inputFile);
            if (inputFileType == InputType::FASTA) {
                *isFastaInput = true;
                PBLOG_WARN << "Input is FASTA. Output BAM file cannot be used for polishing with "
                              "GenomicConsensus!";
            } else if (inputFileType == InputType::FASTQ) {
                *isFastqInput = true;
                PBLOG_WARN << "Input is FASTQ. Output BAM file cannot be used for polishing with "
                              "GenomicConsensus!";
            }
        } else {
            InputType firstFileType = DetermineInputTypeFastx(args[0]);
            InputType secondFileType = DetermineInputTypeFastx(args[1]);
            if (firstFileType == InputType::FASTQ && secondFileType == InputType::FASTQ) {
                PBLOG_FATAL << "Both input files are of type fastq. Please check your inputs.";
                std::exit(EXIT_FAILURE);
            }
            if (firstFileType == InputType::FASTA && secondFileType == InputType::FASTQ) {
                referenceFile = args[0];
                inputFile = args[1];
                *isFastqInput = true;
                PBLOG_WARN << "Input is FASTQ. Output BAM file cannot be used for polishing with "
                              "GenomicConsensus!";
            } else if (firstFileType == InputType::FASTQ && secondFileType == InputType::FASTA) {
                referenceFile = args[1];
                inputFile = args[0];
                *isFastqInput = true;
                PBLOG_WARN << "Input is FASTQ. Output BAM file cannot be used for polishing with "
                              "GenomicConsensus!";
            } else {
                referenceFile = args[0];
                inputFile = args[1];
                *isFastaInput = true;
                PBLOG_WARN << "Input is FASTA. Output BAM file cannot be used for polishing with "
                              "GenomicConsensus!";
            }
        }
    } else if (file0Type == InputType::READ && file1Type == InputType::READ) {
        PBLOG_FATAL << "Both input files are of type READ. Please check your inputs.";
        std::exit(EXIT_FAILURE);
    }
    PBLOG_INFO << "READ input file: " << inputFile;
    PBLOG_INFO << "REF  input file: " << referenceFile;
    if (!Utility::FileExists(inputFile)) {
        PBLOG_FATAL << "Input data file does not exist: " << inputFile;
        std::exit(EXIT_FAILURE);
    }
    *isFromJson = Utility::FileExtension(inputFile) == "json";
    if (*isFromJson) inputFile = UnpackJson(inputFile);
    *isFromXML = Utility::FileExtension(inputFile) == "xml";
    BAM::DataSet dsInput;
    try {
        if (!*isFastaInput && !*isFastqInput) dsInput = BAM::DataSet(inputFile);
    } catch (...) {
        PBLOG_FATAL << UNKNOWN_FILE_TYPES;
        std::exit(EXIT_FAILURE);
    }
    *inputType = dsInput.Type();
    bool fromSubreadset = false;
    bool fromConsensuReadSet = false;
    bool fromTranscriptSet = false;
    const auto IsUnrolled = [&]() {
        bool isUnrolled = settings->AlignMode == AlignmentMode::UNROLLED;
        if (isUnrolled)
            PBLOG_INFO << "Will not automatically set preset based on JSON input, because unrolled "
                          "mode via --zmw or --hqregion has been set!";
        return isUnrolled;
    };
    const auto AlignedInput = [&]() {
        *isAlignedInput = true;
        PBLOG_WARN << "Input is aligned reads. Only primary alignments will be "
                      "respected to allow idempotence!";
    };
    switch (*inputType) {
        case BAM::DataSet::TypeEnum::ALIGNMENT:
            AlignedInput();
        case BAM::DataSet::TypeEnum::SUBREAD: {
            if (*isFromJson && !IsUnrolled()) {
                settings->AlignMode = AlignmentMode::SUBREADS;
                PBLOG_INFO << "Setting to SUBREAD preset";
            }
            fromSubreadset = true;
        } break;
        case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
            AlignedInput();
        case BAM::DataSet::TypeEnum::CONSENSUS_READ: {
            if (*isFromJson && !IsUnrolled()) {
                settings->AlignMode = AlignmentMode::CCS;
                PBLOG_INFO << "Setting to CCS preset";
            }
            fromConsensuReadSet = true;
        } break;
        case BAM::DataSet::TypeEnum::TRANSCRIPT_ALIGNMENT:
            AlignedInput();
        case BAM::DataSet::TypeEnum::TRANSCRIPT: {
            if (*isFromJson && !IsUnrolled()) {
                settings->AlignMode = AlignmentMode::ISOSEQ;
                PBLOG_INFO << "Setting to ISOSEQ preset";
            }
            fromTranscriptSet = true;
        } break;
        case BAM::DataSet::TypeEnum::BARCODE:
        case BAM::DataSet::TypeEnum::REFERENCE:
        default: {
            const auto inType = DetermineInputTypeFastx(inputFile);
            if (inType != InputType::FASTA && inType != InputType::FASTQ) {
                PBLOG_FATAL << "Unsupported input data file " << inputFile << " of type "
                            << BAM::DataSet::TypeToName(dsInput.Type());
                std::exit(EXIT_FAILURE);
            }
        }
    }

    if (!Utility::FileExists(referenceFile)) {
        PBLOG_FATAL << "Input reference file does not exist: " << referenceFile;
        std::exit(EXIT_FAILURE);
    }
    std::string reference;
    if (boost::algorithm::ends_with(referenceFile, ".mmi")) {
        reference = referenceFile;
        PBLOG_INFO
            << "Reference input is an index file. Index parameter override options are disabled!";
    } else {
        BAM::DataSet dsRef(referenceFile);
        switch (dsRef.Type()) {
            case BAM::DataSet::TypeEnum::REFERENCE:
                break;
            case BAM::DataSet::TypeEnum::BARCODE:
            case BAM::DataSet::TypeEnum::SUBREAD:
            case BAM::DataSet::TypeEnum::ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_READ:
            default:
                PBLOG_FATAL << "ERROR: Unsupported reference input file " << referenceFile
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

    if (out == "-" && settings->SplitBySample) {
        PBLOG_FATAL << "Cannot split by sample and use output pipe!";
        std::exit(EXIT_FAILURE);
    }

    if (args.size() == 3) {
        std::string outlc = boost::algorithm::to_lower_copy(out);
        bool isToXML = Utility::FileExtension(outlc) == "xml";

        if (isToXML && (boost::algorithm::ends_with(outlc, ".subreadset.xml") ||
                        boost::algorithm::ends_with(outlc, ".consensusreadset.xml") ||
                        boost::algorithm::ends_with(outlc, ".transcriptset.xml"))) {
            PBLOG_FATAL << "Output has to be an alignment dataset! Please use alignmentset.xml, "
                           "consensusalignmentset.xml, or transcriptalignmentset.xml!";
            std::exit(EXIT_FAILURE);
        }

        bool toAlignmentSet = boost::algorithm::ends_with(outlc, ".alignmentset.xml");
        bool toConsensusAlignmentSet =
            boost::algorithm::ends_with(outlc, ".consensusalignmentset.xml");
        bool toTranscriptAlignmentSet =
            boost::algorithm::ends_with(outlc, ".transcriptalignmentset.xml");

        if (isToXML && !toAlignmentSet && !toConsensusAlignmentSet && !toTranscriptAlignmentSet) {
            PBLOG_FATAL << "Output is XML, but of unknown type! Please use alignmentset.xml, "
                           "consensusalignmentset.xml, or transcriptalignmentset.xml";
            std::exit(EXIT_FAILURE);
        }

        if (*isFromXML && isToXML) {
            std::string outputTypeProvided;
            if (toAlignmentSet)
                outputTypeProvided = "AlignmentSet";
            else if (toConsensusAlignmentSet)
                outputTypeProvided = "ConsensusReadSet";
            else if (toTranscriptAlignmentSet)
                outputTypeProvided = "TranscriptSet";

            if (fromSubreadset && !toAlignmentSet) {
                PBLOG_FATAL << "Unsupported dataset combination! Input SubreadSet with output "
                            << outputTypeProvided
                            << "! Please use AlignmentSet as output XML type!";
                std::exit(EXIT_FAILURE);
            }
            if (fromConsensuReadSet && !toConsensusAlignmentSet) {
                PBLOG_FATAL
                    << "Unsupported dataset combination! Input ConsensusReadSet with output "
                    << outputTypeProvided
                    << "! Please use ConsensusAlignmentSet as output XML type!";
                std::exit(EXIT_FAILURE);
            }
            if (fromTranscriptSet && !toTranscriptAlignmentSet) {
                PBLOG_FATAL << "Unsupported dataset combination! Input TranscriptSet with output "
                            << outputTypeProvided
                            << "! Please use TranscriptAlignmentSet as output XML type!";
                std::exit(EXIT_FAILURE);
            }
        }

        if (isToXML && !*isFromXML)
            PBLOG_WARN << "Input is not a dataset, but output is. Please use dataset input for "
                          "full SMRT Link compatibility!";
    }

    return std::tuple<std::string, std::string, std::string>{inputFile, reference, out};
}  // namespace

std::string CreateDataSet(const BAM::DataSet& dsIn, const std::string& refFile,
                          const bool isFromXML, const std::string& outputFile,
                          const std::string& origOutputFile, std::string* id, size_t numAlignments,
                          size_t numBases)
{
    using BAM::DataSet;
    using DataSetElement = PacBio::BAM::internal::DataSetElement;
    // Input dataset
    const auto GetCollection = [&dsIn](std::string* const name,
                                       std::unique_ptr<DataSetElement>* const collection,
                                       std::string* const tags) {
        *name = dsIn.Name();
        *tags = dsIn.Tags();
        const auto md = dsIn.Metadata();
        if (!md.HasChild("Collections")) return false;
        *collection = std::unique_ptr<DataSetElement>(
            new DataSetElement(std::move(md.Child<DataSetElement>("Collections"))));
        return true;
    };

    std::string datasetName;
    std::string tags;
    std::unique_ptr<DataSetElement> collection;
    const bool hasCollection = GetCollection(&datasetName, &collection, &tags);
    const bool hasName = !datasetName.empty();

    std::string metatype = "PacBio.AlignmentFile.AlignmentBamFile";
    std::string outputType = "alignmentset";
    BAM::DataSet::TypeEnum outputEnum = BAM::DataSet::TypeEnum::ALIGNMENT;

    const auto SetOutputAlignment = [&]() {
        metatype = "PacBio.AlignmentFile.AlignmentBamFile";
        outputType = "alignmentset";
        outputEnum = BAM::DataSet::TypeEnum::ALIGNMENT;
    };

    const auto SetOutputConsensus = [&]() {
        metatype = "PacBio.AlignmentFile.ConsensusAlignmentBamFile";
        outputType = "consensusalignmentset";
        outputEnum = BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT;
    };

    const auto SetOutputTranscript = [&]() {
        metatype = "PacBio.AlignmentFile.TranscriptAlignmentBamFile";
        outputType = "transcriptalignmentset";
        outputEnum = BAM::DataSet::TypeEnum::TRANSCRIPT_ALIGNMENT;
    };

    const auto SetFromDatasetInput = [&]() {
        switch (dsIn.Type()) {
            case BAM::DataSet::TypeEnum::SUBREAD:
                SetOutputAlignment();
                break;
            case BAM::DataSet::TypeEnum::CONSENSUS_READ:
                SetOutputConsensus();
                break;
            case BAM::DataSet::TypeEnum::TRANSCRIPT:
                SetOutputTranscript();
                break;
            default:
                throw std::runtime_error("Unsupported input type");
        }
    };

    if (isFromXML) {
        SetFromDatasetInput();
    } else {
        if (boost::algorithm::ends_with(origOutputFile, ".alignmentset.xml")) {
            SetOutputAlignment();
        } else if (boost::algorithm::ends_with(origOutputFile, ".consensusalignmentset.xml")) {
            SetOutputConsensus();
        } else if (boost::algorithm::ends_with(origOutputFile, ".transcriptalignmentset.xml")) {
            SetOutputTranscript();
        } else if (boost::algorithm::ends_with(origOutputFile, ".json")) {
            SetFromDatasetInput();
        } else {
            PBLOG_FATAL << "Unknown file ending. Please use alignmentset.xml, "
                           "consensusalignmentset.xml, or transcriptalignmentset.xml!";
            std::exit(EXIT_FAILURE);
        }
    }
    DataSet ds(outputEnum);
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
    BAM::FileIndex pbi("PacBio.Index.PacBioIndex", fileName + ".bam.pbi");
    resource.FileIndices().Add(pbi);
    BAM::ExternalResource refResource("PacBio.ReferenceFile.ReferenceFastaFile", refFile);
    resource.ExternalResources().Add(refResource);
    ds.ExternalResources().Add(resource);
    std::string name;
    if (hasName)
        name = datasetName;
    else
        name = fileName;

    name += " (aligned)";
    ds.Name(name);
    ds.TimeStampedName(name + "-" + PacBio::BAM::CurrentTimestamp());

    PacBio::BAM::DataSetMetadata metadata(std::to_string(numAlignments), std::to_string(numBases));
    if (hasCollection) metadata.AddChild(*collection.get());
    ds.Metadata(metadata);

    std::string outputDSFileName = outputFile + "." + outputType + ".xml";
    std::ofstream dsOut(outputDSFileName);
    *id = ds.UniqueId();
    ds.SaveToStream(dsOut);
    return outputDSFileName;
}

std::string OutputFilePrefix(const std::string& outputFile)
{
    // Check if output type is a dataset
    const std::string outputExt = Utility::FileExtension(outputFile);
    std::string prefix = outputFile;

    const std::string outputExtLc = boost::algorithm::to_lower_copy(outputExt);
    if (outputExtLc == "xml") {
        boost::ireplace_last(prefix, ".xml", "");
        boost::ireplace_last(prefix, ".consensusalignmentset", "");
        boost::ireplace_last(prefix, ".alignmentset", "");
        boost::ireplace_last(prefix, ".transcriptalignmentset", "");
    } else if (outputExtLc == "bam") {
        boost::ireplace_last(prefix, ".bam", "");
    } else if (outputExtLc == "json") {
        boost::ireplace_last(prefix, ".json", "");
    } else {
        PBLOG_FATAL << "Unknown file extension for output file: " << outputFile;
        std::exit(EXIT_FAILURE);
    }
    return prefix;
}

static std::string CreateTmpFile(const std::string& outFile)
{
    std::string pipeName;
    int counter = 0;
    do {
        pipeName = std::tmpnam(nullptr);
        if (++counter == 100) break;
        PBLOG_DEBUG << "Trying pipe " << pipeName;
    } while (Utility::FileExists(pipeName));
    if (counter == 100) {
        pipeName = outFile + ".pbmm2.pipe";
        PBLOG_DEBUG << "Trying pipe " << pipeName;
    }
    return pipeName;
}

static void PrintErrorAndAbort(int error)
{
    if (error == EACCES) {
        PBLOG_FATAL << "Pipe error: "
                    << "A component of the path prefix denies search permission, or write "
                       "permission is denied on the parent directory of the FIFO to be "
                       "created.";
    }
    if (error == EEXIST) {
        PBLOG_FATAL << "Pipe error: "
                    << "The named file already exists. Please remove file!";
    }
    if (error == ELOOP) {
        PBLOG_FATAL << "Pipe error: "
                    << "A loop exists in symbolic links encountered during resolution of "
                       "the path argument.";
    }
    if (error == ENAMETOOLONG) {
        PBLOG_FATAL << "Pipe error: "
                    << "The length of the path argument exceeds {PATH_MAX} or a pathname "
                       "component is longer than {NAME_MAX}.";
    }
    if (error == ENOENT) {
        PBLOG_FATAL << "Pipe error: "
                    << "A component of the path prefix specified by path does not name an "
                       "existing directory or path is an empty string.";
    }
    if (error == ENOSPC) {
        PBLOG_FATAL << "Pipe error: "
                    << "The directory that would contain the new file cannot be extended "
                       "or the file system is out of file-allocation resources.";
    }
    if (error == ENOTDIR) {
        PBLOG_FATAL << "Pipe error: "
                    << "A component of the path prefix is not a directory.";
    }
    if (error == EROFS) {
        PBLOG_FATAL << "Pipe error: "
                    << "The named file resides on a read-only file system.";
    }
    if (error == ELOOP) {
        PBLOG_FATAL << "Pipe error: "
                    << "More than {SYMLOOP_MAX} symbolic links were encountered during "
                       "resolution of the path argument.";
    }
    if (error == ENAMETOOLONG) {
        PBLOG_FATAL << "Pipe error: "
                    << "As a result of encountering a symbolic link in resolution of the "
                       "path argument, the length of the substituted pathname string "
                       "exceeded {PATH_MAX}";
    }
    std::exit(EXIT_FAILURE);
}

struct StreamWriter
{
    StreamWriter(BAM::BamHeader header, const std::string& outPrefix, bool sort, int sortThreads,
                 int numThreads, int64_t sortMemory, const std::string& sample = "",
                 const std::string& infix = "")
        : sort_(sort)
        , sample_(sample)
        , outPrefix_(outPrefix)
        , sortThreads_(sortThreads)
        , numThreads_(numThreads)
        , sortMemory_(sortMemory)
        , header_(std::move(header))
    {
        if (outPrefix_ != "-") {
            if (sample_.empty())
                finalOutputPrefix_ = outPrefix_;
            else
                finalOutputPrefix_ = outPrefix_ + '.' + infix;
            finalOutputName_ = finalOutputPrefix_ + ".bam";
        }

        if (!sample_.empty()) {
            auto rgs = header_.ReadGroups();
            header_.ClearReadGroups();
            for (auto& rg : rgs)
                if (rg.Sample() == sample_) header_.AddReadGroup(rg);
        }

        BAM::BamWriter::Config bamWriterConfig;
        std::string outputFile;
        if (sort_) {
            pipeName_ = CreateTmpFile(finalOutputName_);
            int pipe = mkfifo(pipeName_.c_str(), 0666);
            if (pipe == -1) {
                int error = errno;
                PBLOG_FATAL << "Could not open pipe! File name: " << pipeName_;
                PrintErrorAndAbort(error);
            }

            sortThread_ = std::make_unique<std::thread>([&]() {
                int numFiles = 0;
                int numBlocks = 0;
                bam_sort(pipeName_.c_str(), finalOutputName_.c_str(), sortThreads_,
                         sortThreads_ + numThreads_, sortMemory_, &numFiles, &numBlocks);
                PBLOG_INFO << "Merged sorted output from " << numFiles << " files and " << numBlocks
                           << " in-memory blocks";
            });

            outputFile = pipeName_;
            bamWriterConfig.useTempFile = false;
        } else {
            outputFile = finalOutputName_;
            bamWriterConfig.useTempFile = true;
        }
        bamWriter_ = std::make_unique<BAM::BamWriter>(outputFile, header_, bamWriterConfig);
    }

    void Write(const BAM::BamRecord& r) const
    {
        if (!bamWriter_) {
            PBLOG_FATAL << "Nullpointer BamWriter";
            std::exit(EXIT_FAILURE);
        }
        bamWriter_->Write(r);
    }

    int64_t Close()
    {
        bamWriter_.reset();
        if (sort_) {
            Utility::Stopwatch sortTime;
            unlink(pipeName_.c_str());
            if (sortThread_) sortThread_->join();
            return sortTime.ElapsedMilliseconds();
        } else {
            return 0;
        }
    }

    std::string FinalOutputName() { return finalOutputName_; }
    std::string FinalOutputPrefix() { return finalOutputPrefix_; }

private:
    bool sort_;
    std::string sample_;
    std::string outPrefix_;
    int sortThreads_;
    int numThreads_;
    int64_t sortMemory_;
    std::string pipeName_;
    std::unique_ptr<BAM::BamWriter> bamWriter_;
    std::unique_ptr<std::thread> sortThread_;

public:
    std::string finalOutputName_{"-"};
    std::string finalOutputPrefix_{"-"};
    BAM::BamHeader header_;
};
struct StreamWriters
{
    StreamWriters(BAM::BamHeader& header, const std::string& outPrefix, bool splitBySample,
                  bool sort, int sortThreads, int numThreads, int64_t sortMemory)
        : header_(header.DeepCopy())
        , outPrefix_(outPrefix)
        , splitBySample_(splitBySample)
        , sort_(sort)
        , sortThreads_(sortThreads)
        , numThreads_(numThreads)
        , sortMemory_(sortMemory){};

    StreamWriter& at(const std::string& infix, const std::string& sample)
    {
        static const std::string unsplit = "unsplit";
        if (!splitBySample_) {
            if (sampleNameToStreamWriter.find(unsplit) == sampleNameToStreamWriter.cend())
                sampleNameToStreamWriter.emplace(
                    unsplit,
                    std::make_unique<StreamWriter>(header_.DeepCopy(), outPrefix_, sort_,
                                                   sortThreads_, numThreads_, sortMemory_));

            return *sampleNameToStreamWriter.at(unsplit);
        } else {
            if (sampleNameToStreamWriter.find(sample) == sampleNameToStreamWriter.cend())
                sampleNameToStreamWriter.emplace(
                    sample, std::make_unique<StreamWriter>(header_.DeepCopy(), outPrefix_, sort_,
                                                           sortThreads_, numThreads_, sortMemory_,
                                                           sample, infix));
            return *sampleNameToStreamWriter.at(sample);
        }
    }

    std::string WriteDatasetsJson(const BAM::DataSet& inFile, const std::string& origOutFile,
                                  const std::string& refFile, const bool isFromXML,
                                  const bool outputIsJson, const Summary& s,
                                  const std::string& outputFilePrefix, const bool splitSample)
    {
        std::string pbiTiming;
        Timer pbiTimer;
        std::vector<std::string> xmlNames;
        std::vector<std::string> ids;
        for (auto& sample_sw : sampleNameToStreamWriter) {
            BAM::BamFile validationBam(sample_sw.second->FinalOutputName());
            BAM::PbiFile::CreateFrom(validationBam);

            std::string id;
            const auto xmlName =
                CreateDataSet(inFile, refFile, isFromXML, sample_sw.second->FinalOutputPrefix(),
                              origOutFile, &id, s.NumAlns, s.Bases);
            xmlNames.emplace_back(xmlName);
            ids.emplace_back(id);
        }
        pbiTiming = pbiTimer.ElapsedTime();

        if (outputIsJson || splitSample) {
            JSON::Json datastore;
            datastore["createdAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
            datastore["updatedAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
            datastore["version"] = "0.2.2";
            std::vector<JSON::Json> files;
            int i = 0;
            for (auto& sample_sw : sampleNameToStreamWriter) {
                JSON::Json datastoreFile;
                datastoreFile["createdAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
                datastoreFile["description"] = "Aligned and sorted reads as BAM";
                std::ifstream file(sample_sw.second->finalOutputName_,
                                   std::ios::binary | std::ios::ate);
                datastoreFile["fileSize"] = static_cast<int>(file.tellg());

                switch (inFile.Type()) {
                    case BAM::DataSet::TypeEnum::SUBREAD:
                        datastoreFile["fileTypeId"] = "PacBio.DataSet.AlignmentSet";
                        break;
                    case BAM::DataSet::TypeEnum::CONSENSUS_READ:
                        datastoreFile["fileTypeId"] = "PacBio.DataSet.ConsensusAlignmentSet";
                        break;
                    case BAM::DataSet::TypeEnum::TRANSCRIPT:
                        datastoreFile["fileTypeId"] = "PacBio.DataSet.TranscriptAlignmentSet";
                        break;
                    default:
                        throw std::runtime_error("Unsupported input type");
                }

                datastoreFile["isChunked"] = false;
                datastoreFile["modifiedAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
                datastoreFile["name"] = "Aligned reads";
                datastoreFile["path"] = xmlNames[i];
                datastoreFile["sourceId"] = "mapping.tasks.pbmm2_align-out-1";
                boost::uuids::random_generator gen;
                datastoreFile["uniqueId"] = ids[i];
                ++i;
                files.emplace_back(datastoreFile);
            }
            datastore["files"] = files;
            std::ofstream datastoreStream(outputFilePrefix + ".json");
            datastoreStream << datastore.dump(2);
        }
        return pbiTiming;
    }

    std::string Close()
    {
        CreateEmptyIfNoOutput();
        int64_t ms = 0;
        for (auto& sample_sw : sampleNameToStreamWriter) {
            ms += sample_sw.second->Close();
        }
        if (!sort_)
            return "";
        else
            return Timer::ElapsedTimeFromSeconds(ms * 1e6);
    }

    void CreateEmptyIfNoOutput()
    {
        if (sampleNameToStreamWriter.empty()) {
            splitBySample_ = false;
            this->at("", "");
        }
    }

private:
    BAM::BamHeader header_;
    const std::string& outPrefix_;
    bool splitBySample_;
    bool sort_;
    int sortThreads_;
    int numThreads_;
    int64_t sortMemory_;
    std::map<std::string, std::string> movieNameToSampleName;
    std::map<std::string, std::unique_ptr<StreamWriter>> sampleNameToStreamWriter;
};
}  // namespace

int AlignWorkflow::Runner(const CLI::Results& options)
{
    const Timer startTime;
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

    std::string inFileString;
    std::string refFile;
    std::string outFile;
    bool isFromJson = false;
    bool isFromXML = false;
    bool isAlignedInput = false;
    bool isFastaInput = false;
    bool isFastqInput = false;
    BAM::DataSet::TypeEnum inputType;

    std::tie(inFileString, refFile, outFile) =
        CheckPositionalArgs(options.PositionalArguments(), &settings, &isFromJson, &isFromXML,
                            &inputType, &isAlignedInput, &isFastaInput, &isFastqInput);
    BAM::DataSet inFile;
    try {
        if (!isFastaInput && !isFastqInput) inFile = BAM::DataSet(inFileString);
    } catch (...) {
        PBLOG_FATAL << UNKNOWN_FILE_TYPES;
        std::exit(EXIT_FAILURE);
    }

    if (settings.ZMW && !isAlignedInput &&
        (!isFromXML || (inputType != BAM::DataSet::TypeEnum::SUBREAD &&
                        inputType != BAM::DataSet::TypeEnum::ALIGNMENT))) {
        PBLOG_FATAL << "Option --zmw can only be used with a subreadset.xml containing subread + "
                       "scraps BAM files.";
        std::exit(EXIT_FAILURE);
    }
    if (settings.HQRegion && !isAlignedInput &&
        (!isFromXML || (inputType != BAM::DataSet::TypeEnum::SUBREAD &&
                        inputType != BAM::DataSet::TypeEnum::ALIGNMENT))) {
        PBLOG_FATAL
            << "Option --hqregion can only be used with a subreadset.xml containing subread + "
               "scraps BAM files.";
        std::exit(EXIT_FAILURE);
    }

    if (!isFastaInput && !isFastqInput && !settings.Rg.empty()) {
        PBLOG_FATAL << "Cannot override read groups with BAM input. Remove option --rg.";
        std::exit(EXIT_FAILURE);
    }

    std::string outputFilePrefix{"-"};
    bool outputIsXML = false;
    bool outputIsJson = false;
    if (outFile != "-") {
        outputFilePrefix = OutputFilePrefix(outFile);
        const auto ext = Utility::FileExtension(outFile);
        outputIsXML = ext == "xml";
        outputIsJson = ext == "json";
        std::string alnFile = outFile;
        if (outputIsXML || outputIsJson) alnFile = outputFilePrefix + ".bam";

        if (Utility::FileExists(alnFile))
            PBLOG_WARN << "Warning: Overwriting existing output file: " << alnFile;
        if (alnFile != outFile && Utility::FileExists(outFile))
            PBLOG_WARN << "Warning: Overwriting existing output file: " << outFile;
    }

    const FilterFunc filter = [&settings](const AlignedRecord& aln) {
        if (aln.Span <= 0 || aln.Span < settings.MinAlignmentLength) return false;
        if (settings.MinPercConcordance <= 0) return true;
        if (aln.Concordance < settings.MinPercConcordance) return false;
        return true;
    };

    Timer indexTime;
    MM2Helper mm2helper(refFile, settings);
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
        const auto SanitizeSampleName = [](const std::string& in) {
            if (in.empty()) return fallbackSampleName;

            auto trimmed = boost::algorithm::trim_copy(in);
            if (trimmed.empty()) return fallbackSampleName;
            std::string sanitizedName;
            for (const char& c : trimmed) {
                if (c < '!' || c > '~')
                    sanitizedName += '_';
                else
                    sanitizedName += c;
            }
            return sanitizedName;
        };
        const auto SanitizeFileInfix = [](const std::string& in) {
            std::string sanitizedName;
            for (const char& c : in) {
                if (c == '_' || c == '-' || c == '.' || (c >= '0' && c <= '9') ||
                    (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'))
                    sanitizedName += c;
            }
            return sanitizedName;
        };

        std::map<std::string, std::pair<std::string, std::string>> movieNameToSampleAndInfix;
        const auto md = inFile.Metadata();
        if (md.HasChild("Collections")) {
            using DataSetElement = PacBio::BAM::internal::DataSetElement;

            DataSetElement collections = md.Child<DataSetElement>("Collections");
            for (const auto& collectionMetaData : collections.Children()) {
                if (!collectionMetaData.HasAttribute("Context")) {
                    PBLOG_ERROR << "Cannot parse Context attribute of <CollectionMetadata> "
                                   "element. Bailing on biosample parsing.";
                    continue;
                }
                std::string movieName = collectionMetaData.Attribute("Context");
                std::string wellSampleName;
                std::string bioSampleName;
                const auto wellSample = collectionMetaData.Child<DataSetElement>("WellSample");
                if (wellSample.HasAttribute("Name")) wellSampleName = wellSample.Attribute("Name");
                if (wellSample.HasChild("BioSamples")) {
                    const auto bioSamples = wellSample.Child<DataSetElement>("BioSamples");
                    if (bioSamples.HasChild("BioSample")) {
                        const auto bioSample = bioSamples.Child<DataSetElement>("BioSample");
                        if (bioSample.HasAttribute("Name"))
                            bioSampleName = bioSample.Attribute("Name");
                    }
                }
                std::string finalName;
                if (!bioSampleName.empty())
                    finalName = bioSampleName;
                else if (!wellSampleName.empty())
                    finalName = wellSampleName;
                finalName = SanitizeSampleName(finalName);

                movieNameToSampleAndInfix[movieName] = {finalName, SanitizeFileInfix(finalName)};
            }
        }
        std::map<std::string, std::set<std::string>> infixToSamples;
        for (const auto& movie_sampleInfix : movieNameToSampleAndInfix)
            infixToSamples[movie_sampleInfix.second.second].insert(movie_sampleInfix.second.first);

        for (const auto& infix_samples : infixToSamples) {
            auto& infix = infix_samples.first;
            if (infix_samples.second.size() > 1) {
                int infixCounter = 0;
                for (const auto& sample : infix_samples.second) {
                    std::string newInfix = infix + '-' + std::to_string(infixCounter);
                    for (auto& movie_sampleInfix : movieNameToSampleAndInfix) {
                        if (movie_sampleInfix.second.first == sample) {
                            movie_sampleInfix.second.second = newInfix;
                        }
                    }
                    ++infixCounter;
                }
            }
        }

        std::string fastxRgId = "default";
        BAM::BamHeader hdr;
        if ((settings.HQRegion || settings.ZMW) && !isAlignedInput) {
            BAM::ZmwReadStitcher reader(inFile);
            if (reader.HasNext()) {
                auto r = reader.Next();
                hdr = r.Header();
            }
        } else if (!isFastaInput && !isFastqInput) {
            const auto bamFiles = inFile.BamFiles();
            hdr = bamFiles.front().Header();
            for (size_t i = 1; i < bamFiles.size(); ++i)
                hdr += bamFiles.at(i).Header();
        } else {
            std::string rgString = settings.Rg;
            boost::replace_all(rgString, "\\t", "\t");
            if (rgString.empty()) rgString = "@RG\tID:default";
            BAM::ReadGroupInfo rg = BAM::ReadGroupInfo::FromSam(rgString);
            fastxRgId = rg.Id();
            if (rg.MovieName().empty()) rg.MovieName("default");
            hdr.AddReadGroup(rg);
        }
        if (isAlignedInput) {
            hdr.ClearSequences();
        }

        bool performOverrideSampleName = !settings.SampleName.empty();
        std::string overridingSampleName = SanitizeSampleName(settings.SampleName);
        auto rgs = hdr.ReadGroups();
        hdr.ClearReadGroups();
        for (auto& rg : rgs) {
            if (performOverrideSampleName) {
                rg.Sample(overridingSampleName);
            } else if (rg.Sample().empty()) {
                if (isFromXML || isFromJson) {
                    if (movieNameToSampleAndInfix.find(rg.MovieName()) !=
                        movieNameToSampleAndInfix.cend()) {
                        rg.Sample(movieNameToSampleAndInfix[rg.MovieName()].first);
                    } else {
                        PBLOG_INFO << "Cannot find biosample name for movie name " << rg.MovieName()
                                   << "! Will use fallback.";
                        rg.Sample(SanitizeSampleName(""));
                    }
                } else {
                    rg.Sample(SanitizeSampleName(""));
                }
            }
            hdr.AddReadGroup(rg);
        }

        const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
        for (const auto& si : mm2helper.SequenceInfos())
            hdr.AddSequence(si);
        hdr.AddProgram(BAM::ProgramInfo("pbmm2").Name("pbmm2").Version(version).CommandLine(
            "pbmm2 " + settings.CLI));

        PacBio::Parallel::FireAndForget faf(settings.NumThreads, 3);

        writers = std::make_unique<StreamWriters>(hdr, outputFilePrefix, settings.SplitBySample,
                                                  settings.Sort, settings.SortThreads,
                                                  settings.NumThreads, settings.SortMemory);

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
                    writers
                        ->at(movieNameToSampleAndInfix[movieName].second,
                             movieNameToSampleAndInfix[movieName].first)
                        .Write(aln.Record);
                    if (++alignedRecords % settings.ChunkSize == 0) {
                        const auto now = std::chrono::steady_clock::now();
                        auto elapsedSecs =
                            std::chrono::duration_cast<std::chrono::seconds>(now - lastTime)
                                .count();
                        if (elapsedSecs > 5) {
                            lastTime = now;
                            auto elapsedSecTotal =
                                std::chrono::duration_cast<std::chrono::seconds>(now - firstTime)
                                    .count() /
                                60.0;
                            auto alnsPerMin = std::round(alignedReads / elapsedSecTotal);
                            PBLOG_DEBUG << "#Reads, #Aln, #RPM: " << alignedReads << ", "
                                        << alignedRecords << ", " << alnsPerMin;
                        }
                    }
                }
            }
            waiting--;
        };

        const auto FastxToUnalignedBam = [&hdr, &fastxRgId](const std::string& seq,
                                                            const std::string& name,
                                                            const std::string& qual = "") {
            BAM::BamRecord record(hdr);
            record.Impl().SetSequenceAndQualities(seq, qual);
            record.ReadGroupId(fastxRgId);
            record.Impl().Name(name);
            record.Impl().AddTag("qs", 0);
            record.Impl().AddTag("qe", static_cast<int32_t>(seq.size()));
            return record;
        };

        if (isFastaInput) {
            BAM::FastaReader reader(inFileString);
            BAM::FastaSequence fa;
            while (reader.GetNext(fa)) {
                (*records)[i++] = FastxToUnalignedBam(fa.Bases(), fa.Name());
                if (i >= chunkSize) {
                    waiting++;
                    faf.ProduceWith(Submit, std::move(records));
                    records = std::make_unique<std::vector<BAM::BamRecord>>(chunkSize);
                    i = 0;
                }
            }
        } else if (isFastqInput) {
            BAM::FastqReader reader(inFileString);
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
        } else if (isAlignedInput) {
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
    std::string sortTiming = writers->Close();

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
    if (outputIsXML || outputIsJson)
        pbiTiming = writers->WriteDatasetsJson(inFile, outFile, refFile, isFromXML, outputIsJson, s,
                                               outputFilePrefix, settings.SplitBySample);

    PBLOG_INFO << "Mapped Reads: " << alignedReads;
    PBLOG_INFO << "Alignments: " << s.NumAlns;
    PBLOG_INFO << "Mapped Bases: " << s.Bases;
    PBLOG_INFO << "Mean Mapped Concordance: " << meanMappedConcordance << "%";
    PBLOG_INFO << "Max Mapped Read Length: " << maxMappedLength;
    PBLOG_INFO << "Mean Mapped Read Length: " << (1.0 * s.Bases / s.NumAlns);

    PBLOG_INFO << "Index Build/Read Time: " << indexTime.ElapsedTime();
    PBLOG_INFO << "Alignment Time: " << alignmentTime.ElapsedTime();
    if (!sortTiming.empty()) PBLOG_INFO << "Sort Merge Time: " << sortTiming;
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
