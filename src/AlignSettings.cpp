#include "AlignSettings.h"

#include <pbmm2/LibraryInfo.h>
#include "AbortException.h"
#include "Pbmm2GlobalVersion.h"

#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <boost/algorithm/string.hpp>

#include <unistd.h>
#include <cmath>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace PacBio {
namespace minimap2 {
namespace {

int64_t SizeStringToIntMG(const std::string& s)
{
    if (isalpha(s[s.size() - 1])) {
        int64_t size = std::stoll(s.substr(0, s.size() - 1));
        switch (s[s.size() - 1]) {
            case 'm':
            case 'M':
                size <<= 20;
                break;
            case 'g':
            case 'G':
                size <<= 30;
                break;
            default:
                std::ostringstream os;
                os << "Unknown size multiplier " << s[s.size() - 1];
                throw AbortException(os.str());
        }
        return size;
    } else {
        return std::stoll(s);
    }
}

}  // namespace

namespace OptionNames {
// clang-format off

const CLI_v2::Option MinPercConcordance{
R"({
    "names" : ["c", "min-concordance-perc"],
    "description" : "Minimum alignment concordance in percent.",
    "type" : "double",
    "default" : 0,
    "hidden" : true
})"};

const CLI_v2::Option MinPercIdentity{
R"({
    "names" : ["x", "min-id-perc"],
    "description" : "Minimum sequence identity in percent.",
    "type" : "double",
    "default" : 0,
    "hidden" : true
})"};

const CLI_v2::Option MinPercIdentityGapComp{
R"({
    "names" : ["y", "min-gap-comp-id-perc"],
    "description" : "Minimum gap-compressed sequence identity in percent.",
    "type" : "double",
    "default" : 70
})"};

const CLI_v2::Option MinAlignmentLength{
R"({
    "names" : ["l", "min-length"],
    "description" : "Minimum mapped read length in basepairs.",
    "type" : "int",
    "default" : 50
})"};

const CLI_v2::Option SampleName{
R"({
    "names" : ["sample"],
    "description" : [
        "Sample name for all read groups. Defaults, in order of precedence:",
        " SM field in input read group, biosample name, well sample name, \"UnnamedSample\"."
    ],
    "type" : "string"
})"};

const CLI_v2::Option AlignAlignmentModeOpt{
R"({
    "names" : ["preset"],
    "description" : "Set alignment mode. See below for preset parameter details.",
    "type" : "string",
    "choices" : ["SUBREAD", "CCS", "HIFI", "ISOSEQ", "UNROLLED"],
    "default" : "SUBREAD"
})"};

const CLI_v2::Option ChunkSize{
R"({
    "names" : ["chunk-size"],
    "description" : "Process N records per chunk.",
    "type" : "int",
    "default" : 100
})"};

const CLI_v2::Option AlignKmer{
R"({
    "names" : ["k"],
    "description" : "k-mer size (no larger than 28).",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option AlignMinimizerWindowSize{
R"({
    "names" : ["w"],
    "description" : "Minimizer window size.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option GapOpen1{
R"({
    "names" : ["o", "gap-open-1"],
    "description" : "Gap open penalty 1.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option GapOpen2{
R"({
    "names" : ["O", "gap-open-2"],
    "description" : "Gap open penalty 2.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option GapExtension1{
R"({
    "names" : ["e", "gap-extend-1"],
    "description" : "Gap extension penalty 1.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option GapExtension2{
R"({
    "names" : ["E", "gap-extend-2"],
    "description" : "Gap extension penalty 2.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option MatchScore{
R"({
    "names" : ["A"],
    "description" : "Matching score.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option MismatchPenalty{
R"({
    "names" : ["B"],
    "description" : "Mismatch penalty.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option Zdrop{
R"({
    "names" : ["z"],
    "description" : "Z-drop score.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option ZdropInv{
R"({
    "names" : ["Z"],
    "description" : "Z-drop inversion score.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option Bandwidth{
R"({
    "names" : ["r"],
    "description" : "Bandwidth used in chaining and DP-based alignment.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option MaxIntronLength{
R"({
    "names" : ["G"],
    "description" : "Max intron length (changes -r).",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option NonCanon{
R"({
    "names" : ["C"],
    "description" : "Cost for a non-canonical GT-AG splicing (effective in ISOSEQ preset).",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option NoSpliceFlank{
R"({
    "names" : ["no-splice-flank"],
    "description" : "Do not prefer splice flanks GT-AG (effective in ISOSEQ preset)."
})"};

const CLI_v2::Option MedianFilter{
R"({
    "names" : ["median-filter"],
    "description" : "Pick one read per ZMW of median length."
})"};

const CLI_v2::Option Sort{
R"({
    "names" : ["sort"],
    "description" : "Generate sorted BAM file."
})"};

const CLI_v2::Option SortThreadsTC{
R"({
    "names" : ["sort-threads-perc"],
    "description" : [
        "Percentage of threads used exclusively for sorting (absolute number of ",
        "sort threads is capped at 8)."
    ],
    "type" : "int",
    "default" : 25,
    "hidden" : true
})"};

const CLI_v2::Option SortThreads{
R"({
    "names" : ["J", "sort-threads"],
    "description" : "Number of threads used for sorting; 0 means 25% of -j, maximum 8.",
    "type" : "int",
    "default" : 0
})"};

const CLI_v2::Option SortMemory{
R"({
    "names" : ["m", "sort-memory"],
    "description" : "Memory per thread for sorting.",
    "type" : "string",
    "default" : "768M"
})"};

const CLI_v2::Option SortMemoryTC{
R"({
    "names" : ["sort-memory-tc"],
    "description" : "Memory per thread for sorting.",
    "type" : "string",
    "default" : "4G",
    "hidden" : true
})"};

const CLI_v2::Option AlignDisableHPC{
R"({
    "names" : ["u", "no-kmer-compression"],
    "description" : "Disable homopolymer-compressed k-mer (compression is active for SUBREAD & UNROLLED presets)."
})"};

const CLI_v2::Option ZMW{
R"({
    "names" : ["zmw"],
    "description" : "Process ZMW Reads, subreadset.xml input required (activates UNROLLED preset)."
})"};

const CLI_v2::Option HQRegion{
R"({
    "names" : ["hqregion"],
    "description" : "Process HQ region of each ZMW, subreadset.xml input required (activates UNROLLED preset)."
})"};

const CLI_v2::Option Strip{
R"({
    "names" : ["strip"],
    "description" : "Remove all kinetic and extra QV tags. Output cannot be polished."
})"};

const CLI_v2::Option SplitBySample{
R"({
    "names" : ["split-by-sample"],
    "description" : "One output BAM per sample."
})"};

const CLI_v2::Option Rg{
R"({
    "names" : ["rg"],
    "description" : "Read group header line such as '@RG\\tID:xyz\\tSM:abc'. Only for FASTA/Q inputs.",
    "type" : "string"
})"};

const CLI_v2::Option CreatePbi{
R"({
    "names" : ["pbi"],
    "description" : "Generate PBI for BAM only output.",
    "hidden" : true
})"};

const CLI_v2::Option LongJoinFlankRatio{
R"({
    "names" : ["L", "lj-min-ratio"],
    "description" : "Long join flank ratio.",
    "type" : "float",
    "default" : -1
})"};

const CLI_v2::Option NoBAI{
R"({
    "names" : ["no-bai"],
    "description" : "Omit BAI generation for sorted output.",
    "hidden" : true
})"};

const CLI_v2::Option NoTrimming{
R"({
    "names" : ["no-rmt"],
    "description" : "Disable repeated matches trimming.",
    "hidden" : true
})"};

const CLI_v2::Option OutputUnmapped{
R"({
    "names" : ["unmapped"],
    "description" : "Include unmapped records in output."
})"};

const CLI_v2::Option MaxNumAlns{
R"({
    "names" : ["N", "best-n"],
    "description" : "Output at maximum N alignments for each read, 0 means no maximum.",
    "type" : "int",
    "default" : 0
})"};

const CLI_v2::Option CompressSequenceHomopolymers{
R"({
    "names" : ["collapse-homopolymers"],
    "description" : "Collapse homopolymers in reads and reference."
})"};

const CLI_v2::Option MaxGap{
R"({
    "names" : ["g"],
    "description" : "Stop chain enlongation if there are no minimizers in N bp.",
    "type" : "int",
    "default" : -1
})"};

const CLI_v2::Option BamIndexInput{
R"({
    "names" : ["bam-index"],
    "description" : "Generate index for sorted BAM output.",
    "type" : "string",
    "choices" : ["NONE", "BAI", "CSI"],
    "default" : "BAI"
})"};

const CLI_v2::Option EnforcedMapping{
R"({
    "names" : ["enforced-mapping"],
    "description" : "Enforce mapping of reads to reference.",
    "type" : "string",
    "hidden" : true
})"};

const CLI_v2::Option MaxSecondaryAlns{
R"({
    "names" : ["max-secondary-alns"],
    "description" : "Retain at most N secondary alignments prior filtering.",
    "type" : "int",
    "default" : 5,
    "hidden" : true
})"};

const CLI_v2::Option ReportFileJson{
R"({
    "names" : ["report-json"],
    "description" : [
        "Write a mapping stats report JSON to this filename."
    ],
    "type" : "string",
    "hidden": true
})"};

const CLI_v2::Option ShortSACigar{
R"({
    "names" : ["short-sa-cigar"],
    "description" : "Populate SA tag with short cigar representation.",
    "type" : "bool"
})"};


const CLI_v2::PositionalArgument Reference {
R"({
    "name" : "ref.fa|xml|mmi",
    "description" : "Reference FASTA, ReferenceSet XML, or Reference Index"
})"};

const CLI_v2::PositionalArgument Input {
R"({
    "name" : "in.bam|xml|fa|fq|gz|fofn",
    "description" : "Input BAM, DataSet XML, FASTA, or FASTQ"
})"};

const CLI_v2::PositionalArgument Output {
R"({
    "name" : "out.aligned.bam|xml",
    "description" : "Output BAM or DataSet XML",
    "required" : false
})"};

// clang-format on
}  // namespace OptionNames

AlignSettings::AlignSettings(const PacBio::CLI_v2::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , MinPercConcordance(options[OptionNames::MinPercConcordance])
    , MinPercIdentity(options[OptionNames::MinPercIdentity])
    , MinPercIdentityGapComp(options[OptionNames::MinPercIdentityGapComp])
    , MinAlignmentLength(options[OptionNames::MinAlignmentLength])
    , SampleName(options[OptionNames::SampleName])
    , ReportFileJson(options[OptionNames::ReportFileJson])
    , ChunkSize(options[OptionNames::ChunkSize])
    , MedianFilter(options[OptionNames::MedianFilter])
    , Sort(options[OptionNames::Sort])
    , ZMW(options[OptionNames::ZMW])
    , HQRegion(options[OptionNames::HQRegion])
    , Strip(options[OptionNames::Strip])
    , SplitBySample(options[OptionNames::SplitBySample])
    , Rg(options[OptionNames::Rg])
    , CreatePbi(options[OptionNames::CreatePbi])
    , OutputUnmapped(options[OptionNames::OutputUnmapped])
    , CompressSequenceHomopolymers(options[OptionNames::CompressSequenceHomopolymers])
{
    MM2Settings::Kmer = options[OptionNames::AlignKmer];
    MM2Settings::MinimizerWindowSize = options[OptionNames::AlignMinimizerWindowSize];
    MM2Settings::GapOpen1 = options[OptionNames::GapOpen1];
    MM2Settings::GapOpen2 = options[OptionNames::GapOpen2];
    MM2Settings::GapExtension1 = options[OptionNames::GapExtension1];
    MM2Settings::GapExtension2 = options[OptionNames::GapExtension2];
    MM2Settings::MatchScore = options[OptionNames::MatchScore];
    MM2Settings::MismatchPenalty = options[OptionNames::MismatchPenalty];
    MM2Settings::Zdrop = options[OptionNames::Zdrop];
    MM2Settings::ZdropInv = options[OptionNames::ZdropInv];
    MM2Settings::Bandwidth = options[OptionNames::Bandwidth];
    MM2Settings::MaxIntronLength = options[OptionNames::MaxIntronLength];
    MM2Settings::NonCanon = options[OptionNames::NonCanon];
    MM2Settings::NoSpliceFlank = options[OptionNames::NoSpliceFlank];
    MM2Settings::DisableHPC = options[OptionNames::AlignDisableHPC];
    MM2Settings::LongJoinFlankRatio = options[OptionNames::LongJoinFlankRatio];
    MM2Settings::NoTrimming = options[OptionNames::NoTrimming];
    MM2Settings::MaxNumAlns = options[OptionNames::MaxNumAlns];
    MM2Settings::MaxGap = options[OptionNames::MaxGap];
    MM2Settings::EnforcedMapping = std::string(options[OptionNames::EnforcedMapping]);
    if (!MM2Settings::EnforcedMapping.empty()) MM2Settings::NoTrimming = true;
    MM2Settings::MaxSecondaryAlns = options[OptionNames::MaxSecondaryAlns];
    MM2Settings::ShortSACigar = options[OptionNames::ShortSACigar];

    const bool noBai = options[OptionNames::NoBAI];
    const std::string bamIdx = options[OptionNames::BamIndexInput];

    int numAvailableCores = std::thread::hardware_concurrency();
    const unsigned int rawRequestedNThreads = options[PacBio::CLI_v2::Builtin::NumThreads];
    int requestedNThreads =
        static_cast<int>(rawRequestedNThreads);  // since we're doing subtractions here
    std::string requestedMemory = options[OptionNames::SortMemory];
    int sortThreadPerc = 25;
    SortThreads = options[OptionNames::SortThreads];

    if (sortThreadPerc > 50)
        PBLOG_WARN << "Please allocate less than 50% of threads for sorting. Currently allocated: "
                   << sortThreadPerc << "%!";

    if (requestedNThreads > numAvailableCores) {
        PBLOG_WARN << "Requested more threads for alignment (" << requestedNThreads
                   << ") than system-wide available (" << numAvailableCores << ")";
    }

    int availableThreads = ThreadCount(requestedNThreads);
    if (Sort) {
        if (SortThreads == 0) {
            SortThreads = std::min(
                std::max(static_cast<int>(std::round(availableThreads * sortThreadPerc / 100.0)),
                         1),
                8);
            MM2Settings::NumThreads = std::max(availableThreads - SortThreads, 1);
        } else if (SortThreads != 0) {
            if (requestedNThreads == 0)
                availableThreads = std::max(availableThreads - SortThreads, 1);
            if (availableThreads + SortThreads > numAvailableCores) {
                PBLOG_WARN << "Requested more threads for sorting (" << SortThreads
                           << ") and alignment (" << availableThreads
                           << ") than system-wide available (" << numAvailableCores << ")";
            }
            MM2Settings::NumThreads = availableThreads;
            if (SortThreads > numAvailableCores) {
                PBLOG_WARN << "Requested more threads for sorting (" << SortThreads
                           << ") than system-wide available (" << numAvailableCores << ")!";
                SortThreads = ThreadCount(SortThreads);
            }
            while (MM2Settings::NumThreads + SortThreads > numAvailableCores ||
                   (MM2Settings::NumThreads == 1 && SortThreads == 1)) {
                SortThreads = std::max(SortThreads - 1, 1);
                if (MM2Settings::NumThreads + SortThreads <= numAvailableCores ||
                    (MM2Settings::NumThreads == 1 && SortThreads == 1))
                    break;
                MM2Settings::NumThreads = std::max(MM2Settings::NumThreads - 1, 1);
            }
        }
    } else {
        MM2Settings::NumThreads = availableThreads;
    }
    SortMemory = SizeStringToIntMG(requestedMemory);

    if (!Sort) {
        if (SortThreads != 0)
            PBLOG_WARN
                << "Requested " << SortThreads
                << " threads for sorting, without specifying --sort. Please check your input.";
        const std::string pureMemory = options[OptionNames::SortMemory];
        if (pureMemory != "768M")
            PBLOG_WARN
                << "Requested " << pureMemory
                << " memory for sorting, without specifying --sort. Please check your input.";
    }

    if (Sort) {
        std::string suffix;
        const auto MemoryToHumanReadable = [](int64_t memInBytes, float* roundedMemory,
                                              std::string* suffix) {
            static constexpr int64_t ninek = 1000;
            *roundedMemory = memInBytes;
            if (memInBytes >> 10 < ninek) {
                *roundedMemory /= 1000;
                *suffix = "K";
            } else if (memInBytes >> 20 < ninek) {
                *roundedMemory = (memInBytes >> 10) / 1024.0;
                *suffix = "M";
            } else if (memInBytes >> 30 < ninek) {
                *roundedMemory = (memInBytes >> 20) / 1024.0;
                *suffix = "G";
            }
        };
        float maxMemSortFloat;
        std::string maxMemSortSuffix;
        MemoryToHumanReadable(SortMemory * SortThreads, &maxMemSortFloat, &maxMemSortSuffix);
        PBLOG_INFO << "Using " << MM2Settings::NumThreads << " threads for alignments, "
                   << SortThreads << " threads for sorting, and " << maxMemSortFloat
                   << maxMemSortSuffix << " bytes RAM for sorting.";

        BamIdx = BamIndex::_from_string(bamIdx.c_str());

        if (noBai) {
            PBLOG_WARN << "Overriding --bam-index with --no-bai!";
            BamIdx = BamIndex::NONE;
        }
    } else {
        PBLOG_INFO << "Using " << MM2Settings::NumThreads << " threads for alignments.";
    }

    const std::map<std::string, AlignmentMode> alignModeMap{{"SUBREAD", AlignmentMode::SUBREADS},
                                                            {"ISOSEQ", AlignmentMode::ISOSEQ},
                                                            {"CCS", AlignmentMode::CCS},
                                                            {"HIFI", AlignmentMode::CCS},
                                                            {"UNROLLED", AlignmentMode::UNROLLED}};

    const std::string alignModeUsr = options[OptionNames::AlignAlignmentModeOpt];
    const std::string alingModeUpr = boost::to_upper_copy(alignModeUsr);
    if (alignModeMap.find(alingModeUpr) == alignModeMap.cend()) {
        throw AbortException("Could not find --preset " + alignModeUsr);
    }
    MM2Settings::AlignMode = alignModeMap.at(alingModeUpr);
    int inputFilterCounts = ZMW + MedianFilter + HQRegion;
    if (inputFilterCounts > 1) {
        throw AbortException(
            "Options --zmw, --hqregion and --median-filter are mutually exclusive.");
    }
    if (ZMW || HQRegion) {
        if (ChunkSize != 100)
            PBLOG_WARN << "Cannot change --chunk-size in --zmw/--hqregion mode. Parameters "
                          "--chunk-size is forced to 1.";
        ChunkSize = 1;
        MM2Settings::AlignMode = AlignmentMode::UNROLLED;
    }

    if (!Rg.empty() && !boost::contains(Rg, "ID") && !boost::starts_with(Rg, "@RG\t")) {
        throw AbortException(
            "Invalid @RG line. Missing ID field. Please provide following format: "
            "'@RG\\tID:xyz\\tSM:abc'");
    }
    if (MM2Settings::LongJoinFlankRatio > 1) {
        throw AbortException("Option -L,--lj-min-ratio has to be between a ratio betweem 0 and 1.");
    }

    if (!Sort && noBai) {
        PBLOG_WARN << "Option --no-bai has no effect without option --sort!";
    }

    if (MM2Settings::GapOpen1 < -1 || MM2Settings::GapOpen2 < -1 ||
        MM2Settings::GapExtension1 < -1 || MM2Settings::GapExtension2 < -1) {
        throw AbortException("Gap options have to be strictly positive.");
    }
    if (MM2Settings::Kmer < -1 || MM2Settings::Kmer == 0 || MM2Settings::MinimizerWindowSize < -1 ||
        MM2Settings::MinimizerWindowSize == 0) {
        throw AbortException("Index parameter -k and -w must be positive.");
    }

    if (MM2Settings::MaxNumAlns < 0) {
        throw AbortException("Parameter --best-n, -N must be positive.");
    }

    // Override Sample Name for all Read Groups, disable SplitBySample.
    if (!SampleName.empty() && SplitBySample) {
        PBLOG_WARN << "Options --split-by-sample and --sample are mutually exclusive. Option "
                      "--sample will be applied and --split-by-sample is ignored!";
        SplitBySample = false;
    }
}

int32_t AlignSettings::ThreadCount(int32_t n)
{
    const int32_t m = std::thread::hardware_concurrency();
    if (n <= 0) n = m + n;  // permit n <= 0 to subtract from max threads
    return std::max(1, std::min(m, n));
}

PacBio::CLI_v2::Interface AlignSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pbmm2 align", "Align PacBio reads to reference sequences",
                                Pbmm2::LibraryInfo().Release};

    i.Example("pbmm2 align ref.referenceset.xml movie.subreadset.xml ref.movie.alignmentset.xml");

    // clang-format off
    i.AddPositionalArguments({
        OptionNames::Reference,
        OptionNames::Input,
        OptionNames::Output
    });

    i.AddOptionGroup("Basic Options", {
        OptionNames::ChunkSize,
        OptionNames::NoTrimming,

        // hidden
        OptionNames::SortMemoryTC,
        OptionNames::CreatePbi,
        OptionNames::EnforcedMapping,
        OptionNames::MaxSecondaryAlns,
    });

    i.AddOptionGroup("Sorting Options", {
        OptionNames::Sort,
        OptionNames::SortMemory,
        OptionNames::SortThreads,
        OptionNames::SortThreadsTC,
    });

    i.AddOptionGroup("Parameter Set Options", {
        OptionNames::AlignAlignmentModeOpt,
    });

    i.AddOptionGroup("General Parameter Override Options", {
        OptionNames::AlignKmer,
        OptionNames::AlignMinimizerWindowSize,
        OptionNames::AlignDisableHPC,
        OptionNames::MatchScore,
        OptionNames::MismatchPenalty,
        OptionNames::Zdrop,
        OptionNames::ZdropInv,
        OptionNames::Bandwidth,
        OptionNames::MaxGap,
    });

    i.AddOptionGroup("Gap Parameter Override Options (a k-long gap costs min{o+k*e,O+k*E})", {
        OptionNames::GapOpen1,
        OptionNames::GapOpen2,
        OptionNames::GapExtension1,
        OptionNames::GapExtension2,
        OptionNames::LongJoinFlankRatio,
    });

    i.AddOptionGroup("IsoSeq Parameter Override Options", {
        OptionNames::MaxIntronLength,
        OptionNames::NonCanon,
        OptionNames::NoSpliceFlank,
    });

    i.AddOptionGroup("Read Group Options", {
        OptionNames::SampleName,
        OptionNames::Rg,
    });

    i.AddOptionGroup("Identity Filter Options (combined with AND)", {
        OptionNames::MinPercConcordance,
        OptionNames::MinPercIdentity,
        OptionNames::MinPercIdentityGapComp,
    });

    i.AddOptionGroup("Output Options", {
        OptionNames::MinAlignmentLength,
        OptionNames::MaxNumAlns,
        OptionNames::Strip,
        OptionNames::SplitBySample,
        OptionNames::OutputUnmapped,
        OptionNames::BamIndexInput,
        OptionNames::NoBAI,
        OptionNames::ShortSACigar,
        OptionNames::ReportFileJson
    });

    i.AddOptionGroup("Input Manipulation Options (mutually exclusive)", {
        OptionNames::MedianFilter,
        OptionNames::ZMW,
        OptionNames::HQRegion,
    });

    i.AddOptionGroup("Sequence Manipulation Options", {
        OptionNames::CompressSequenceHomopolymers
    });

    i.HelpFooter(R"(Alignment modes of --preset:
    SUBREAD     : -k 19 -w 10    -o 5 -O 56 -e 4 -E 1 -A 2 -B 5 -z 400 -Z 50  -r 2000   -L 0.5 -g 5000
    CCS or HiFi : -k 19 -w 10 -u -o 5 -O 56 -e 4 -E 1 -A 2 -B 5 -z 400 -Z 50  -r 2000   -L 0.5 -g 5000
    ISOSEQ      : -k 15 -w 5  -u -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 200000 -L 0.5 -g 2000 -C 5 -G 200000
    UNROLLED    : -k 15 -w 15    -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 2000   -L 0.5 -g 10000)");

    i.RegisterVersionPrinter(Pbmm2::PrintPbmm2VersionSingle);

    // clang-format on
    return i;
}
}  // namespace minimap2
}  // namespace PacBio
