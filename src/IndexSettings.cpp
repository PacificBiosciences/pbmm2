// Author: Armin Töpfer
#include "IndexSettings.h"

#include <map>

#include <boost/algorithm/string.hpp>

#include <pbmm2/Pbmm2Version.h>

#include "AbortException.h"

namespace PacBio {
namespace minimap2 {
namespace OptionNames {
// clang-format off

const CLI_v2::Option IndexAlignmentModeOpt{
R"({
    "names" : ["preset"],
    "description" : "Set alignment mode. See below for preset parameter details.",
    "type" : "string",
    "choices" : ["SUBREAD", "CCS", "ISOSEQ", "UNROLLED"],
    "default" : "SUBREAD"
})"};

const CLI_v2::Option IndexKmer{
R"({
    "names" : ["k"],
    "description" : "k-mer size (no larger than 28).",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option IndexMinimizerWindowSize{
R"({
    "names" : ["w"],
    "description" : "Minimizer window size.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option IndexDisableHPC{
R"({
    "names" : ["u", "no-kmer-compression"],
    "description" : "Disable homopolymer-compressed k-mer (compression is active for SUBREAD & UNROLLED presets)."
})"};

const CLI_v2::PositionalArgument InputReference {
R"({
    "name" : "ref.fa|xml",
    "description" : "Reference FASTA, ReferenceSet XML"
})"};

const CLI_v2::PositionalArgument OutputIndex {
R"({
    "name" : "out.mmi",
    "description" : "Output Reference Index"
})"};

// clang-format on
}  // namespace OptionNames

IndexSettings::IndexSettings(const PacBio::CLI_v2::Results& options)
    : CLI{options.InputCommandLine()}, InputFiles{options.PositionalArguments()}
{
    MM2Settings::Kmer = options[OptionNames::IndexKmer];
    MM2Settings::MinimizerWindowSize = options[OptionNames::IndexMinimizerWindowSize];
    MM2Settings::DisableHPC = options[OptionNames::IndexDisableHPC];
    MM2Settings::NumThreads = options.NumThreads();

    const std::map<std::string, AlignmentMode> alignModeMap{{"SUBREAD", AlignmentMode::SUBREADS},
                                                            {"ISOSEQ", AlignmentMode::ISOSEQ},
                                                            {"CCS", AlignmentMode::CCS},
                                                            {"HIFI", AlignmentMode::CCS},
                                                            {"UNROLLED", AlignmentMode::UNROLLED}};

    const std::string alignModeUsr = options[OptionNames::IndexAlignmentModeOpt];
    const std::string alingModeUpr = boost::to_upper_copy(alignModeUsr);
    if (alignModeMap.find(alingModeUpr) == alignModeMap.cend()) {
        throw AbortException("Could not find --preset " + alignModeUsr);
    }
    MM2Settings::AlignMode = alignModeMap.at(alingModeUpr);

    if (MM2Settings::Kmer < -1 || MM2Settings::Kmer == 0 || MM2Settings::MinimizerWindowSize < -1 ||
        MM2Settings::MinimizerWindowSize == 0) {
        throw AbortException("Index parameter -k and -w must be positive.");
    }
}

PacBio::CLI_v2::Interface IndexSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pbmm2 index", "Index reference and store as .mmi file",
                                PacBio::Pbmm2FormattedVersion()};

    i.Example("pbmm2 index ref.referenceset.xml ref.mmi");

    // clang-format off
    i.AddPositionalArguments({
        OptionNames::InputReference,
        OptionNames::OutputIndex
    });

    i.AddOptionGroup("Parameter Set Option", {
        OptionNames::IndexAlignmentModeOpt,
    });

    i.AddOptionGroup("Parameter Override Options", {
        OptionNames::IndexKmer,
        OptionNames::IndexMinimizerWindowSize,
        OptionNames::IndexDisableHPC
    });

    i.HelpFooter(R"(Alignment modes of --preset:
    SUBREAD     : -k 19 -w 10
    CCS or HiFi : -k 19 -w 10 -u
    ISOSEQ      : -k 15 -w 5  -u
    UNROLLED    : -k 15 -w 15
    )");

    // clang-format on
    return i;
}
}  // namespace minimap2
}  // namespace PacBio
