// Author: Armin Töpfer

#include <fstream>

#include <pbbam/DataSet.h>
#include <pbcopper/utility/FileUtils.h>
#include <boost/algorithm/string.hpp>

#include "AlignSettings.h"

#include "InputOutputUX.h"

namespace PacBio {
namespace minimap2 {
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

InputType InputOutputUX::DetermineInputTypeFastx(const std::string& in)
{
    char firstChar = ' ';
    std::ifstream f(in);
    f.read(&firstChar, 1);
    if (firstChar == '>')
        return InputType::FASTA;
    else if (firstChar == '@')
        return InputType::FASTQ;
    else if (firstChar == 'M')
        return InputType::MMI;
    else {
        PBLOG_FATAL << "Unkown file type for file " << in << " starting with character "
                    << firstChar;
        std::exit(EXIT_FAILURE);
    }
}

InputType InputOutputUX::DetermineInputTypeApprox(std::string inputFile)
{
    if (!Utility::FileExists(inputFile)) {
        PBLOG_FATAL << "Input data file does not exist: " << inputFile;
        std::exit(EXIT_FAILURE);
    }
    auto fileExt = Utility::FileExtension(inputFile);
    if (fileExt == "json") {
        inputFile = UnpackJson(inputFile);
        fileExt = Utility::FileExtension(inputFile);
    }
    if (fileExt == "mmi") return InputType::MMI;
    if (fileExt == "fq") return InputType::FASTX;
    if (fileExt == "fastq") return InputType::FASTX;
    if (fileExt == "fa") return InputType::FASTX;
    if (fileExt == "fasta") return InputType::FASTX;
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
            return InputType::BAM;
        case BAM::DataSet::TypeEnum::BARCODE:
        case BAM::DataSet::TypeEnum::REFERENCE:
            return InputType::FASTX;
        default:
            PBLOG_FATAL << "Unsupported input data file " << inputFile << " of type "
                        << BAM::DataSet::TypeToName(dsInput.Type());
            std::exit(EXIT_FAILURE);
    }
}

UserIO InputOutputUX::CheckPositionalArgs(const std::vector<std::string>& args,
                                          AlignSettings& settings)
{
    UserIO uio;
    if (args.size() < 2) {
        PBLOG_FATAL << "Please provide at least the input arguments: reference input output!";
        PBLOG_FATAL << "EXAMPLE: pbmm2 reference.fasta input.subreads.bam output.bam";
        std::exit(EXIT_FAILURE);
    }

    std::string inputFile;
    std::string referenceFile;
    const auto file0Type = InputOutputUX::DetermineInputTypeApprox(args[0]);
    const auto file1Type = InputOutputUX::DetermineInputTypeApprox(args[1]);
    const bool firstRefFastx = file0Type == InputType::FASTX || file0Type == InputType::MMI;
    const bool secondRefFastx = file1Type == InputType::FASTX || file1Type == InputType::MMI;

    if (file0Type == InputType::BAM && secondRefFastx) {  // BAM <FASTX|MMI>
        inputFile = args[0];
        referenceFile = args[1];
    } else if (firstRefFastx && file1Type == InputType::BAM) {  // <FASTX|MMI> BAM
        referenceFile = args[0];
        inputFile = args[1];
    } else if (firstRefFastx && secondRefFastx) {  // <FASTX|MMI> <FASTX|MMI>
        const bool firstMmi = file0Type == InputType::MMI;
        const bool secondMmi = file1Type == InputType::MMI;
        const int32_t numMmi = firstMmi + secondMmi;
        const std::string noGC =
            "Output BAM file cannot be used for polishing with GenomicConsensus!";

        if (numMmi == 2) {  // MMI MMI
            PBLOG_FATAL << "Both input files are of type MMI. Please check your inputs.";
            std::exit(EXIT_FAILURE);
        } else if (numMmi == 1) {  // <FASTX|MMI> <FASTX|MMI>
            if (firstMmi) {        // MMI FASTX
                referenceFile = args[0];
                inputFile = args[1];
            } else {  //  FASTX MMI
                referenceFile = args[1];
                inputFile = args[0];
            }
            InputType inputFileType = InputOutputUX::DetermineInputTypeFastx(inputFile);
            if (inputFileType == InputType::FASTA) {  // FASTX == FASTA
                uio.isFastaInput = true;
                PBLOG_WARN << "Input is FASTA." << noGC;
            } else if (inputFileType == InputType::FASTQ) {  // FASTX == FASTQ
                uio.isFastqInput = true;
                PBLOG_WARN << "Input is FASTQ." << noGC;
            }
        } else {  // FASTX FASTX
            InputType firstFileType = DetermineInputTypeFastx(args[0]);
            InputType secondFileType = DetermineInputTypeFastx(args[1]);
            const bool firstFastq = firstFileType == InputType::FASTQ;
            const bool secondFastq = secondFileType == InputType::FASTQ;
            const bool firstFasta = firstFileType == InputType::FASTA;
            const bool secondFasta = secondFileType == InputType::FASTA;
            if (firstFastq && secondFastq) {  // FASTQ FASTQ
                PBLOG_FATAL << "Both input files are of type FASTQ. Please check your inputs.";
                std::exit(EXIT_FAILURE);
            }
            if (firstFasta && secondFastq) {  // FASTA FASTQ
                referenceFile = args[0];
                inputFile = args[1];
                uio.isFastqInput = true;
                PBLOG_WARN << "Input is FASTQ." << noGC;
            } else if (firstFastq && secondFasta) {  // FASTQ FASTA
                referenceFile = args[1];
                inputFile = args[0];
                uio.isFastqInput = true;
                PBLOG_WARN << "Input is FASTQ." << noGC;
            } else {  // FASTA FASTA
                referenceFile = args[0];
                inputFile = args[1];
                uio.isFastaInput = true;
                PBLOG_WARN << "Input is FASTA." << noGC;
            }
        }
    } else if (file0Type == InputType::BAM && file1Type == InputType::BAM) {  // BAM BAM
        PBLOG_FATAL << "Both input files are of type BAM. Please check your inputs.";
        std::exit(EXIT_FAILURE);
    }
    PBLOG_INFO << "READ input file: " << inputFile;
    PBLOG_INFO << "REF  input file: " << referenceFile;

    if (InputOutputUX::DetermineInputTypeFastx(referenceFile) == InputType::FASTQ) {
        PBLOG_FATAL << "Cannot use FASTQ input as reference. Please use FASTA!";
        std::exit(EXIT_FAILURE);
    }

    auto inputFileExt = Utility::FileExtension(inputFile);
    if (inputFileExt == "json") {
        uio.isFromJson = true;
        inputFile = UnpackJson(inputFile);
        inputFileExt = Utility::FileExtension(inputFile);
    }
    uio.isFromXML = inputFileExt == "xml";

    if (!uio.isFastaInput && !uio.isFastqInput) {
        BAM::DataSet dsInput = BAM::DataSet(inputFile);
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
                [[fallthrough]];
            case BAM::DataSet::TypeEnum::SUBREAD: {
                if (uio.isFromJson && !IsUnrolled()) {
                    settings.AlignMode = AlignmentMode::SUBREADS;
                    PBLOG_INFO << "Setting to SUBREAD preset";
                }
                uio.isFromSubreadset = true;
            } break;
            case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
                AlignedInput();
                [[fallthrough]];
            case BAM::DataSet::TypeEnum::CONSENSUS_READ: {
                if (uio.isFromJson && !IsUnrolled()) {
                    settings.AlignMode = AlignmentMode::CCS;
                    PBLOG_INFO << "Setting to CCS preset";
                }
                uio.isFromConsensuReadSet = true;
            } break;
            case BAM::DataSet::TypeEnum::TRANSCRIPT_ALIGNMENT:
                AlignedInput();
                [[fallthrough]];
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
                const auto inType = DetermineInputTypeFastx(inputFile);
                if (inType != InputType::FASTA && inType != InputType::FASTQ) {
                    PBLOG_FATAL << "Unsupported input data file " << inputFile << " of type "
                                << BAM::DataSet::TypeToName(dsInput.Type());
                    std::exit(EXIT_FAILURE);
                }
            }
        }
    }

    std::string reference;
    if (Utility::FileExtension(referenceFile) == "mmi") {
        reference = referenceFile;
        PBLOG_INFO << "Reference input is an index file. Index parameter override options are "
                      "disabled!";
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

    if (args.size() == 3)
        uio.outFile = args[2];
    else
        uio.outFile = "-";

    if (uio.outFile == "-" && settings.SplitBySample) {
        PBLOG_FATAL << "Cannot split by sample and use output pipe!";
        std::exit(EXIT_FAILURE);
    }

    if (args.size() == 3) {
        const std::string outlc = boost::algorithm::to_lower_copy(uio.outFile);
        const auto outExt = Utility::FileExtension(outlc);
        uio.isToXML = outExt == "xml";
        uio.isToJson = outExt == "json";

        if (uio.isToXML && (boost::algorithm::ends_with(outlc, ".subreadset.xml") ||
                            boost::algorithm::ends_with(outlc, ".consensusreadset.xml") ||
                            boost::algorithm::ends_with(outlc, ".transcriptset.xml"))) {
            PBLOG_FATAL << "Output has to be an alignment dataset! Please use alignmentset.xml, "
                           "consensusalignmentset.xml, or transcriptalignmentset.xml!";
            std::exit(EXIT_FAILURE);
        }

        const bool toAlignmentSet = boost::algorithm::ends_with(outlc, ".alignmentset.xml");
        const bool toConsensusAlignmentSet =
            boost::algorithm::ends_with(outlc, ".consensusalignmentset.xml");
        const bool toTranscriptAlignmentSet =
            boost::algorithm::ends_with(outlc, ".transcriptalignmentset.xml");

        if (uio.isToXML && !toAlignmentSet && !toConsensusAlignmentSet &&
            !toTranscriptAlignmentSet) {
            PBLOG_FATAL << "Output is XML, but of unknown type! Please use alignmentset.xml, "
                           "consensusalignmentset.xml, or transcriptalignmentset.xml";
            std::exit(EXIT_FAILURE);
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
                PBLOG_FATAL << "Unsupported dataset combination! Input SubreadSet with output "
                            << outputTypeProvided
                            << "! Please use AlignmentSet as output XML type!";
                std::exit(EXIT_FAILURE);
            }
            if (uio.isFromConsensuReadSet && !toConsensusAlignmentSet) {
                PBLOG_FATAL
                    << "Unsupported dataset combination! Input ConsensusReadSet with output "
                    << outputTypeProvided
                    << "! Please use ConsensusAlignmentSet as output XML type!";
                std::exit(EXIT_FAILURE);
            }
            if (uio.isFromTranscriptSet && !toTranscriptAlignmentSet) {
                PBLOG_FATAL << "Unsupported dataset combination! Input TranscriptSet with output "
                            << outputTypeProvided
                            << "! Please use TranscriptAlignmentSet as output XML type!";
                std::exit(EXIT_FAILURE);
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
    } else {
        uio.outPrefix = '-';
    }

    uio.inFile = inputFile;
    uio.refFile = reference;
    return uio;
}

std::string InputOutputUX::CreateDataSet(const BAM::DataSet& dsIn, const std::string& refFile,
                                         const bool isFromXML, const std::string& outputFile,
                                         const std::string& origOutputFile, std::string* id,
                                         size_t numAlignments, size_t numBases)
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
        PBLOG_FATAL << "Unknown file extension for output file: " << outputFile;
        std::exit(EXIT_FAILURE);
    }
    return prefix;
}
}  // namespace minimap2
}  // namespace PacBio
