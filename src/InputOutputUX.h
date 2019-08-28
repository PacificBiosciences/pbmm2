// Author: Armin Töpfer

#pragma once

#include <string>
#include <vector>

#include <pbcopper/json/JSON.h>

#include <pbbam/DataSet.h>

namespace PacBio {
namespace minimap2 {
namespace {
static const std::string UNKNOWN_FILE_TYPES =
    "Could not determine read input type(s). Please do not mix data types, such as BAM+FASTQ. File "
    "of files may only contain BAMs or datasets.";
}

struct AlignSettings;

enum class InputType : int
{
    BAM = 0,
    XML,
    XML_BAM,
    XML_FASTA,
    FASTX,
    MMI,
    FASTA,
    FASTQ,
    DATASET,
    FOFN_BAM,
    FOFN_FASTA,
    FOFN_FASTQ,
    UNKNOWN
};

struct UserIO
{
    std::string inFile;
    std::string refFile;
    std::string outFile;
    std::string outPrefix;
    std::string unpackedFromJson;
    bool isFromFofn = false;
    bool isFromJson = false;
    bool isFromXML = false;
    bool isToXML = false;
    bool isToJson = false;
    bool isAlignedInput = false;
    bool isFastaInput = false;
    bool isFastqInput = false;
    bool isFromSubreadset = false;
    bool isFromConsensuReadSet = false;
    bool isFromTranscriptSet = false;
    bool isFromMmi = false;
    BAM::DataSet::TypeEnum inputType;
    std::vector<std::string> inputFiles;
};

class InputOutputUX
{
public:
    // Dumps json file content into a Json object
    static JSON::Json ReadJson(const std::string& jsonInputFile);

    // Takes JSON path as input and returns the file path of the datastore
    static std::string UnpackJson(const std::string& jsonInputFile);

    // Guess file type based on first character for fasta and fastq
    static InputType DetermineInputTypeFastx(const std::string& in);

    // Based on input file, guess file type on file suffix or
    // dataset type.
    static InputType DetermineInputTypeApprox(std::string inputFile, UserIO& uio);

    static UserIO CheckPositionalArgs(const std::vector<std::string>& args,
                                      AlignSettings& settings);

    static std::string CreateDataSet(const BAM::DataSet& dsIn, const std::string& refFile,
                                     const bool isFromXML, const std::string& outputFile,
                                     const std::string& origOutputFile, std::string* id,
                                     size_t numAlignments, size_t numBases);

    static std::string OutPrefix(const std::string& outputFile);
};
}  // namespace minimap2
}  // namespace PacBio
