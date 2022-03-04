#include "InputOutputUX.h"

#include "AbortException.h"
#include "AlignSettings.h"

#include <pbbam/DataSet.h>
#include <pbcopper/utility/FileUtils.h>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <sstream>

namespace PacBio {
namespace minimap2 {
namespace {
bool InputTypeEquality(const InputType& t0, const InputType& t1)
{
    if (t0 == t1) return true;
    if (t0 == InputType::BAM && t1 == InputType::XML_BAM) return true;
    if (t1 == InputType::BAM && t0 == InputType::XML_BAM) return true;
    return false;
}
InputType DetermineInputFileSuffix(const std::string& inputFile)
{
    using TypeEnum = BAM::DataSet::TypeEnum;
    if (boost::iends_with(inputFile, "fq") || boost::iends_with(inputFile, "fastq") ||
        boost::iends_with(inputFile, "fq.gz") || boost::iends_with(inputFile, "fastq.gz"))
        return InputType::FASTQ;

    if (boost::iends_with(inputFile, "fa") || boost::iends_with(inputFile, "fasta") ||
        boost::iends_with(inputFile, "fa.gz") || boost::iends_with(inputFile, "fasta.gz") ||
        boost::iends_with(inputFile, "fsa"))
        return InputType::FASTA;

    if (boost::iends_with(inputFile, "bam")) return InputType::BAM;

    if (boost::iends_with(inputFile, "xml")) {
        BAM::DataSet dsInput;
        try {
            dsInput = BAM::DataSet{inputFile};
        } catch (...) {
            throw AbortException(UNKNOWN_FILE_TYPES);
        }
        switch (dsInput.Type()) {
            case TypeEnum::ALIGNMENT:
            case TypeEnum::SUBREAD:
            case TypeEnum::CONSENSUS_ALIGNMENT:
            case TypeEnum::CONSENSUS_READ:
            case TypeEnum::TRANSCRIPT_ALIGNMENT:
            case TypeEnum::TRANSCRIPT:
                return InputType::XML_BAM;
                break;
            case TypeEnum::BARCODE:
            case TypeEnum::REFERENCE:
                return InputType::XML_FASTA;
                break;
            default: {
                std::ostringstream os;
                os << "Unsupported input data file " << inputFile << " of type "
                   << BAM::DataSet::TypeToName(dsInput.Type());
                throw AbortException(os.str());
            }
        }
    }

    if (boost::iends_with(inputFile, "mmi")) return InputType::MMI;

    throw AbortException("Unknown file suffix of " + inputFile);
}
InputType DetermineFofnContent(const std::string& fofnInputFile, UserIO& uio)
{
    std::ifstream infile(fofnInputFile);
    std::unique_ptr<InputType> type;
    std::string line;
    while (std::getline(infile, line)) {
        boost::trim(line);
        const InputType t = DetermineInputFileSuffix(line);
        if (!Utility::FileExists(line)) {
            throw AbortException("Input fofn contains non-existing file: " + line);
        }
        if (!type) {
            type = std::make_unique<InputType>(t);
        } else if (!InputTypeEquality(*type, t)) {
            throw AbortException(
                "Input fofn contains different file types. This is not supported.");
        }
        uio.inputFiles.emplace_back(line);
    }
    switch (*type) {
        case InputType::BAM:
        case InputType::XML_BAM:
            return InputType::FOFN_BAM;
        case InputType::FASTA:
            return InputType::FOFN_FASTA;
        case InputType::FASTQ:
            return InputType::FOFN_FASTQ;
        default:
            throw AbortException("Unsupported file types in file " + fofnInputFile);
    }
}
}  // namespace
JSON::Json InputOutputUX::ReadJson(const std::string& jsonInputFile)
{
    std::ifstream ifs(jsonInputFile);
    JSON::Json j;
    ifs >> j;
    return j;
}
std::string InputOutputUX::UnpackJson(const std::string& jsonInputFile)
{
    JSON::Json j = ReadJson(jsonInputFile);
    std::string inputFile;
    const auto panic = [](const std::string& error) {
        throw AbortException("JSON Datastore: " + error);
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

InputType InputOutputUX::DetermineInputTypeApprox(std::string inputFile, UserIO& uio)
{
    if (!Utility::FileExists(inputFile)) {
        throw AbortException("Input data file does not exist: " + inputFile);
    }
    if (boost::iends_with(inputFile, "json")) inputFile = UnpackJson(inputFile);
    if (boost::iends_with(inputFile, "fofn")) return DetermineFofnContent(inputFile, uio);

    return DetermineInputFileSuffix(inputFile);
}

UserIO InputOutputUX::CheckPositionalArgs(const std::vector<std::string>& args,
                                          AlignSettings& settings)
{
    UserIO uio;
    if (args.size() < 2) {
        throw AbortException(
            "Please provide at least the input arguments: reference input output!\n"
            "EXAMPLE: pbmm2 reference.fasta input.subreads.bam output.bam");
    }

    std::string inputFile;
    std::string referenceFile;
    const auto file0Type = InputOutputUX::DetermineInputTypeApprox(args[0], uio);
    const auto file1Type = InputOutputUX::DetermineInputTypeApprox(args[1], uio);

    const auto IsExactCombination = [&](const InputType& type0, const InputType& type1) {
        return (file0Type == type0 && file1Type == type1);
    };

    if (file0Type == InputType::MMI || file1Type == InputType::MMI) uio.isFromMmi = true;

    const std::string noGC = " Output BAM file cannot be used for polishing with GenomicConsensus!";
    if (IsExactCombination(InputType::BAM, InputType::BAM)) {
        throw AbortException("Both input files are of type BAM. Please check your inputs.");
    } else if (IsExactCombination(InputType::FASTQ, InputType::FASTQ)) {
        throw AbortException("Both input files are of type FASTQ. Please check your inputs.");
    } else if (IsExactCombination(InputType::MMI, InputType::MMI)) {
        throw AbortException("Both input files are of type MMI. Please check your inputs.");
    } else if (IsExactCombination(InputType::XML_BAM, InputType::XML_BAM)) {
        throw AbortException(
            "Both input files are of type BAM from XML. Please check your inputs.");
    } else if (IsExactCombination(InputType::XML_FASTA, InputType::XML_FASTA)) {
        throw AbortException(
            "Both input files are of type FASTA from XML. Please check your inputs.");
    } else if (IsExactCombination(InputType::FOFN_BAM, InputType::FOFN_BAM)) {
        throw AbortException(
            "Both input files are of type BAM from FOFN. Please check your inputs.");
    } else if (IsExactCombination(InputType::FOFN_FASTA, InputType::FOFN_FASTA)) {
        throw AbortException(
            "Both input files are of type FASTA from FOFN. Please check your inputs.");
    } else if (IsExactCombination(InputType::FOFN_FASTQ, InputType::FOFN_FASTQ)) {
        throw AbortException(
            "Both input files are of type FASTQ from FOFN. Please check your inputs.");
    } else if (IsExactCombination(InputType::FASTA, InputType::FASTA) ||
               IsExactCombination(InputType::XML_FASTA, InputType::FASTA) ||
               IsExactCombination(InputType::MMI, InputType::FASTA)) {
        PBLOG_WARN << "Input is FASTA." << noGC;
        uio.isFastaInput = true;
        referenceFile = args[0];
        inputFile = args[1];
    } else if (IsExactCombination(InputType::FASTA, InputType::FOFN_FASTA) ||
               IsExactCombination(InputType::XML_FASTA, InputType::FOFN_FASTA) ||
               IsExactCombination(InputType::MMI, InputType::FOFN_FASTA)) {
        PBLOG_WARN << "Input is FASTA FOFN." << noGC;
        uio.isFromFofn = true;
        uio.isFastaInput = true;
        referenceFile = args[0];
        inputFile = args[1];
    } else if (IsExactCombination(InputType::MMI, InputType::XML_BAM) ||
               IsExactCombination(InputType::MMI, InputType::BAM) ||
               IsExactCombination(InputType::FASTA, InputType::XML_BAM) ||
               IsExactCombination(InputType::FASTA, InputType::BAM) ||
               IsExactCombination(InputType::XML_FASTA, InputType::XML_BAM) ||
               IsExactCombination(InputType::XML_FASTA, InputType::BAM)) {
        referenceFile = args[0];
        inputFile = args[1];
    } else if (IsExactCombination(InputType::MMI, InputType::FOFN_BAM) ||
               IsExactCombination(InputType::FASTA, InputType::FOFN_BAM) ||
               IsExactCombination(InputType::XML_FASTA, InputType::FOFN_BAM)) {
        uio.isFromFofn = true;
        referenceFile = args[0];
        inputFile = args[1];
    } else if (IsExactCombination(InputType::MMI, InputType::FASTQ) ||
               IsExactCombination(InputType::FASTA, InputType::FASTQ) ||
               IsExactCombination(InputType::XML_FASTA, InputType::FASTQ)) {
        PBLOG_WARN << "Input is FASTQ." << noGC;
        uio.isFastqInput = true;
        referenceFile = args[0];
        inputFile = args[1];
    } else if (IsExactCombination(InputType::MMI, InputType::FOFN_FASTQ) ||
               IsExactCombination(InputType::FASTA, InputType::FOFN_FASTQ) ||
               IsExactCombination(InputType::XML_FASTA, InputType::FOFN_FASTQ)) {
        PBLOG_WARN << "Input is FASTQ FOFN." << noGC;
        uio.isFromFofn = true;
        uio.isFastqInput = true;
        referenceFile = args[0];
        inputFile = args[1];
    } else if (IsExactCombination(InputType::XML_BAM, InputType::MMI) ||
               IsExactCombination(InputType::BAM, InputType::MMI) ||
               IsExactCombination(InputType::XML_BAM, InputType::FASTA) ||
               IsExactCombination(InputType::BAM, InputType::FASTA) ||
               IsExactCombination(InputType::XML_BAM, InputType::XML_FASTA) ||
               IsExactCombination(InputType::BAM, InputType::XML_FASTA)) {
        inputFile = args[0];
        referenceFile = args[1];
    } else if (IsExactCombination(InputType::FOFN_BAM, InputType::MMI) ||
               IsExactCombination(InputType::FOFN_BAM, InputType::FASTA) ||
               IsExactCombination(InputType::FOFN_BAM, InputType::XML_FASTA)) {
        uio.isFromFofn = true;
        inputFile = args[0];
        referenceFile = args[1];
    } else if (IsExactCombination(InputType::FASTQ, InputType::MMI) ||
               IsExactCombination(InputType::FASTQ, InputType::XML_FASTA) ||
               IsExactCombination(InputType::FASTQ, InputType::FASTA)) {
        PBLOG_WARN << "Input is FASTQ." << noGC;
        uio.isFastqInput = true;
        inputFile = args[0];
        referenceFile = args[1];
    } else if (IsExactCombination(InputType::FOFN_FASTQ, InputType::MMI) ||
               IsExactCombination(InputType::FOFN_FASTQ, InputType::XML_FASTA) ||
               IsExactCombination(InputType::FOFN_FASTQ, InputType::FASTA)) {
        PBLOG_WARN << "Input is FASTQ FOFN." << noGC;
        uio.isFastqInput = true;
        uio.isFromFofn = true;
        inputFile = args[0];
        referenceFile = args[1];
    } else if (IsExactCombination(InputType::FASTA, InputType::MMI) ||
               IsExactCombination(InputType::FASTA, InputType::XML_FASTA)) {
        PBLOG_WARN << "Input is FASTA." << noGC;
        uio.isFastaInput = true;
        inputFile = args[0];
        referenceFile = args[1];
    } else if (IsExactCombination(InputType::FOFN_FASTA, InputType::MMI) ||
               IsExactCombination(InputType::FOFN_FASTA, InputType::XML_FASTA)) {
        PBLOG_WARN << "Input is FASTA FOFN." << noGC;
        uio.isFastaInput = true;
        uio.isFromFofn = true;
        inputFile = args[0];
        referenceFile = args[1];
    } else {
        throw AbortException("Unknown combination");
    }

    if (uio.inputFiles.empty()) {
        uio.inputFiles.emplace_back(inputFile);
    }

    PBLOG_INFO << "READ input file: " << inputFile;
    PBLOG_INFO << "REF  input file: " << referenceFile;

    auto inputFileExt = Utility::FileExtension(inputFile);
    if (inputFileExt == "json") {
        uio.isFromJson = true;
        inputFile = UnpackJson(inputFile);
        uio.unpackedFromJson = inputFile;
        inputFileExt = Utility::FileExtension(inputFile);
    }
    uio.isFromXML = inputFileExt == "xml";

    if (!uio.isFastaInput && !uio.isFastqInput) {
        BAM::DataSet dsInput;
        try {
            dsInput = BAM::DataSet{inputFile};
        } catch (...) {
            throw AbortException(UNKNOWN_FILE_TYPES);
        }
        uio.inputType = dsInput.Type();

        const auto IsUnrolled = [&]() {
            bool isUnrolled = settings.AlignMode == AlignmentMode::UNROLLED;
            if (isUnrolled)
                PBLOG_INFO << "Will not automatically set preset based on JSON input, because "
                              "unrolled "
                              "mode via --zmw or --hqregion has been set!";
            return isUnrolled;
        };
        const auto AlignedInput = [&]() {
            uio.isAlignedInput = true;
            PBLOG_WARN << "Input is aligned reads. Only primary alignments will be "
                          "respected to allow idempotence!";
        };
        switch (uio.inputType) {
            case BAM::DataSet::TypeEnum::ALIGNMENT:
                AlignedInput();
                /* Falls through. */
            case BAM::DataSet::TypeEnum::SUBREAD: {
                if (uio.isFromJson && !IsUnrolled()) {
                    settings.AlignMode = AlignmentMode::SUBREADS;
                    PBLOG_INFO << "Setting to SUBREAD preset";
                }
                uio.isFromSubreadset = true;
            } break;
            case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
                AlignedInput();
                /* Falls through. */
            case BAM::DataSet::TypeEnum::CONSENSUS_READ: {
                if (uio.isFromJson && !IsUnrolled()) {
                    settings.AlignMode = AlignmentMode::CCS;
                    PBLOG_INFO << "Setting to CCS preset";
                }
                uio.isFromConsensuReadSet = true;
            } break;
            case BAM::DataSet::TypeEnum::TRANSCRIPT_ALIGNMENT:
                AlignedInput();
                /* Falls through. */
            case BAM::DataSet::TypeEnum::TRANSCRIPT: {
                if (uio.isFromJson && !IsUnrolled()) {
                    settings.AlignMode = AlignmentMode::ISOSEQ;
                    PBLOG_INFO << "Setting to ISOSEQ preset";
                }
                uio.isFromTranscriptSet = true;
            } break;
            case BAM::DataSet::TypeEnum::BARCODE:
            case BAM::DataSet::TypeEnum::REFERENCE:
            default: {
                std::ostringstream os;
                os << "Unsupported input data file " << inputFile << " of type "
                   << BAM::DataSet::TypeToName(dsInput.Type());
                throw AbortException(os.str());
            }
        }
    }

    std::string reference;
    if (Utility::FileExtension(referenceFile) == "mmi") {
        reference = referenceFile;
        PBLOG_INFO << "Reference input is an index file. Index parameter override options are "
                      "disabled!";
    } else {
        BAM::DataSet dsRef;
        try {
            dsRef = BAM::DataSet(referenceFile);
        } catch (...) {
            throw AbortException(UNKNOWN_FILE_TYPES);
        }
        if (dsRef.Type() != BAM::DataSet::TypeEnum::REFERENCE) {
            std::ostringstream os;
            os << "ERROR: Unsupported reference input file " << referenceFile << " of type "
               << BAM::DataSet::TypeToName(dsRef.Type());
            throw AbortException(os.str());
        }
        const auto fastaFiles = dsRef.FastaFiles();
        if (fastaFiles.size() != 1) {
            std::ostringstream os;
            os << "Only one reference sequence allowed!";
            throw AbortException(os.str());
        }
        reference = fastaFiles.front();
    }

    if (DetermineInputFileSuffix(reference) == InputType::FASTQ) {
        throw AbortException("Cannot use FASTQ input as reference. Please use FASTA!");
    }

    if (args.size() == 3)
        uio.outFile = args[2];
    else
        uio.outFile = "-";

    if (uio.outFile == "-" && settings.SplitBySample) {
        throw AbortException("Cannot split by sample and use output pipe!");
    }

    if (uio.outFile == "-" && settings.CreatePbi) {
        throw AbortException("Cannot generate pbi and use output pipe!");
    }

    if (args.size() == 3) {
        const std::string outlc = boost::algorithm::to_lower_copy(uio.outFile);
        const auto outExt = Utility::FileExtension(outlc);
        uio.isToXML = outExt == "xml";
        uio.isToJson = outExt == "json";

        if (uio.isToXML && (boost::algorithm::ends_with(outlc, ".subreadset.xml") ||
                            boost::algorithm::ends_with(outlc, ".consensusreadset.xml") ||
                            boost::algorithm::ends_with(outlc, ".transcriptset.xml"))) {
            throw AbortException(
                "Output has to be an alignment dataset! Please use alignmentset.xml, "
                "consensusalignmentset.xml, or transcriptalignmentset.xml!");
        }

        const bool toAlignmentSet = boost::algorithm::ends_with(outlc, ".alignmentset.xml");
        const bool toConsensusAlignmentSet =
            boost::algorithm::ends_with(outlc, ".consensusalignmentset.xml");
        const bool toTranscriptAlignmentSet =
            boost::algorithm::ends_with(outlc, ".transcriptalignmentset.xml");

        if (uio.isToXML && !toAlignmentSet && !toConsensusAlignmentSet &&
            !toTranscriptAlignmentSet) {
            throw AbortException(
                "Output is XML, but of unknown type! Please use alignmentset.xml, "
                "consensusalignmentset.xml, or transcriptalignmentset.xml");
        }

        if (uio.isFromXML && uio.isToXML) {
            std::string outputTypeProvided;
            if (toAlignmentSet)
                outputTypeProvided = "AlignmentSet";
            else if (toConsensusAlignmentSet)
                outputTypeProvided = "ConsensusReadSet";
            else if (toTranscriptAlignmentSet)
                outputTypeProvided = "TranscriptSet";

            if (uio.isFromSubreadset && !toAlignmentSet) {
                throw AbortException(
                    "Unsupported dataset combination! Input SubreadSet with output " +
                    outputTypeProvided + "! Please use AlignmentSet as output XML type!");
            }
            if (uio.isFromConsensuReadSet && !toConsensusAlignmentSet) {
                throw AbortException(
                    "Unsupported dataset combination! Input ConsensusReadSet with output " +
                    outputTypeProvided + "! Please use ConsensusAlignmentSet as output XML type!");
            }
            if (uio.isFromTranscriptSet && !toTranscriptAlignmentSet) {
                throw AbortException(
                    "Unsupported dataset combination! Input TranscriptSet with output " +
                    outputTypeProvided + "! Please use TranscriptAlignmentSet as output XML type!");
            }
        }

        if (uio.isToXML && !uio.isFromXML)
            PBLOG_WARN << "Input is not a dataset, but output is. Please use dataset input for "
                          "full SMRT Link compatibility!";

        uio.outPrefix = OutPrefix(uio.outFile);
        std::string alnFile = uio.outFile;
        if (uio.isToXML || uio.isToJson) alnFile = uio.outPrefix + ".bam";

        if (Utility::FileExists(alnFile))
            PBLOG_WARN << "Warning: Overwriting existing output file: " << alnFile;
        if (alnFile != uio.outFile && Utility::FileExists(uio.outFile))
            PBLOG_WARN << "Warning: Overwriting existing output file: " << uio.outFile;
    } else if (args.size() == 2) {
        uio.outPrefix = '-';
    } else if (args.size() == 1 || args.size() > 3) {
        throw AbortException("Incorrect number of arguments. Accepted are at most three!");
    }

    uio.inFile = inputFile;
    uio.refFile = reference;

    if (settings.Sort && settings.BamIdx != BamIndex::_from_index(0)) {
        uio.bamIndex = settings.BamIdx;
    }
    return uio;
}

std::string InputOutputUX::CreateDataSet(const BAM::DataSet& dsIn, const std::string& refFile,
                                         const bool isFromXML, const std::string& outputFile,
                                         const std::string& origOutputFile, std::string* id,
                                         size_t numAlignments, size_t numBases,
                                         const BamIndex& bamIndex)
{
    using BAM::DataSet;
    using TypeEnum = BAM::DataSet::TypeEnum;
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
    TypeEnum outputEnum = TypeEnum::ALIGNMENT;

    const auto SetOutputAlignment = [&]() {
        metatype = "PacBio.AlignmentFile.AlignmentBamFile";
        outputType = "alignmentset";
        outputEnum = TypeEnum::ALIGNMENT;
    };

    const auto SetOutputConsensus = [&]() {
        metatype = "PacBio.AlignmentFile.ConsensusAlignmentBamFile";
        outputType = "consensusalignmentset";
        outputEnum = TypeEnum::CONSENSUS_ALIGNMENT;
    };

    const auto SetOutputTranscript = [&]() {
        metatype = "PacBio.AlignmentFile.TranscriptAlignmentBamFile";
        outputType = "transcriptalignmentset";
        outputEnum = TypeEnum::TRANSCRIPT_ALIGNMENT;
    };

    const auto SetFromDatasetInput = [&]() {
        switch (dsIn.Type()) {
            case TypeEnum::SUBREAD:
                SetOutputAlignment();
                break;
            case TypeEnum::CONSENSUS_READ:
                SetOutputConsensus();
                break;
            case TypeEnum::TRANSCRIPT:
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
            throw AbortException(
                "Unknown file ending. Please use alignmentset.xml, "
                "consensusalignmentset.xml, or transcriptalignmentset.xml!");
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
    BAM::ExternalResource resource(metatype, outputFile + ".bam");
    BAM::FileIndex pbi("PacBio.Index.PacBioIndex", outputFile + ".bam.pbi");
    resource.FileIndices().Add(pbi);

    switch (bamIndex) {
        case BamIndex::BAI: {
            BAM::FileIndex index("PacBio.Index.BamIndex", outputFile + ".bam.bai");
            resource.FileIndices().Add(index);
            break;
        }
        case BamIndex::CSI: {
            BAM::FileIndex index("PacBio.Index.CsiIndex", outputFile + ".bam.csi");
            resource.FileIndices().Add(index);
            break;
        }
        default:
            break;
    }

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

std::string InputOutputUX::OutPrefix(const std::string& outputFile)
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
        throw AbortException("Unknown file extension for output file: " + outputFile);
    }
    return prefix;
}
}  // namespace minimap2
}  // namespace PacBio
