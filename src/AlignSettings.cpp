// Author: Armin Töpfer

#include <map>

#include <Pbmm2Version.h>

#include "AlignSettings.h"

namespace PacBio {
namespace minimap2 {
namespace OptionNames {
// clang-format off
static const CLI::Option HelpOption{
    "help",
    {"h","help"},
    "Output this help.",
    CLI::Option::BoolType()
};
static const CLI::Option VersionOption{
    "version",
    {"version"},
    "Output version information.",
    CLI::Option::BoolType()
};
static const CLI::Option LogLevelOption{
    "log_level",
    {"log-level"},
    R"(Set log level: "TRACE", "DEBUG", "INFO", "WARN", "FATAL".)",
    CLI::Option::StringType("WARN"),
    {"TRACE", "DEBUG", "INFO", "WARN", "FATAL"}
};
const PlainOption LogFile{
    "log_file",
    { "log-file" },
    "Log to a File",
    "Log to a file, instead of stdout.",
    CLI::Option::StringType("")
};
const PlainOption NumThreads{
    "numthreads",
    { "j", "num-threads" },
    "Number of Threads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption MinAccuracy{
    "minaccuracy",
    { "min-accuracy" },
    "Minimum Alignment Accuracy",
    "Minimum alignment accuracy (not effective in ISOSEQ preset).",
    CLI::Option::FloatType(0.75f)
};
const PlainOption MinAlignmentLength{
    "minalnlength",
    { "min-length" },
    "Minimum Alignment Length",
    "Minimum alignment length.",
    CLI::Option::IntType(50)
};
const PlainOption SampleName{
    "biosample_name",
    { "sample-name" },
    "Sample Name",
    "Add sample name to read groups.",
    CLI::Option::StringType()
};
const PlainOption AlignModeOpt{
    "align_mode",
    { "preset" },
    "Alignment mode",
    "Set alignment mode:\n  - \"SUBREAD\" -k 19 -w 10 -d 5 -i 56 -D 4 -I 1 -A 2 -B 5 -z 400 -Z 50 -r 2000\n"
    "  - \"ISOSEQ\" -k 15 -w 5 -d 2 -i 32 -D 1 -I 0 -A 1 -B 2 -z 200 -Z 100 -C 5 -r 200000 -G 200000\n"
    "Default",
    CLI::Option::StringType("SUBREAD"),
    {"SUBREAD", "ISOSEQ"}
};
const PlainOption ChunkSize{
    "chunk_size",
    { "chunk-size" },
    "Chunk Size",
    "Process N records per chunk.",
    CLI::Option::IntType(100)
};
const PlainOption Kmer{
    "kmer_size",
    { "k" },
    "K-mer Size",
    "k-mer size (no larger than 28).",
    CLI::Option::IntType(-1)
};
const PlainOption MinimizerWindowSize{
    "minimizer_window_size",
    { "w" },
    "Minizer Window Size",
    "Minizer window size.",
    CLI::Option::IntType(-1)
};
const PlainOption GapOpenDelete{
    "gap_open_del",
    { "d" },
    "Deletion Gap Open Penalty",
    "Deletion gap open penalty.",
    CLI::Option::IntType(-1)
};
const PlainOption GapOpenInsert{
    "gap_open_ins",
    { "i" },
    "Insertion Gap Open Penalty",
    "Insertion gap open penalty.",
    CLI::Option::IntType(-1)
};
const PlainOption GapExtensionDelete{
    "gap_extension_del",
    { "D" },
    "Deletion Gap Extension Penalty",
    "Deletion gap extension penalty.",
    CLI::Option::IntType(-1)
};
const PlainOption GapExtensionInsert{
    "gap_extension_ins",
    { "I" },
    "Insertion Gap Extension Penalty",
    "Insertion gap extension penalty.",
    CLI::Option::IntType(-1)
};
const PlainOption MatchScore{
    "match_score",
    { "A" },
    "Matching Score",
    "Matching score.",
    CLI::Option::IntType(-1)
};
const PlainOption MismatchPenalty{
    "mismatch_penalty",
    { "B" },
    "Mismatch Penalty",
    "Mismatch penalty.",
    CLI::Option::IntType(-1)
};
const PlainOption Zdrop{
    "zdrop_score",
    { "z" },
    "Z-Drop Score",
    "Z-drop score.",
    CLI::Option::IntType(-1)
};
const PlainOption ZdropInv{
    "zdrop_score_ins",
    { "Z" },
    "Z-Drop Inversion Score",
    "Z-drop inversion score.",
    CLI::Option::IntType(-1)
};
const PlainOption Bandwidth{
    "bandwidth",
    { "r" },
    "Bandwidth Used in Chaining and DP-based Alignment",
    "Bandwidth used in chaining and DP-based alignment.",
    CLI::Option::IntType(-1)
};
const PlainOption MaxIntronLength{
    "max_intron_length",
    { "G" },
    "Max Intron Length (effective in ISOSEQ mode; changing bandwidth)",
    "Max intron length (changes -r).",
    CLI::Option::IntType(-1)
};
const PlainOption NonCanon{
    "non_canon",
    { "C" },
    "Cost For a Non-Canonical GT-AG Splicing (effective in ISOSEQ mode)",
    "Cost for a non-canonical GT-AG splicing.",
    CLI::Option::IntType(-1)
};
const PlainOption NoSpliceFlank{
    "no_splice_flank",
    { "no-splice-flank" },
    "Do Not Prefer Splice Flanks GT-AG (effective in ISOSEQ mode)",
    "Do not prefer splice flanks GT-AG.",
    CLI::Option::BoolType(false)
};
const PlainOption MedianFilter{
    "median_filter",
    { "median-filter" },
    "Pick One Read per ZMW of Median Length",
    "Pick one read per ZMW of median length.",
    CLI::Option::BoolType(false)
};
// clang-format on
}  // namespace OptionNames

AlignSettings::AlignSettings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , MinAccuracy(options[OptionNames::MinAccuracy])
    , MinAlignmentLength(options[OptionNames::MinAlignmentLength])
    , LogFile{options[OptionNames::LogFile].get<decltype(LogFile)>()}
    , LogLevel{options.LogLevel()}
    , SampleName{options[OptionNames::SampleName].get<decltype(SampleName)>()}
    , ChunkSize(options[OptionNames::ChunkSize])
    , MedianFilter(options[OptionNames::MedianFilter])
{
    MM2Settings::Kmer = options[OptionNames::Kmer];
    MM2Settings::MinimizerWindowSize = options[OptionNames::MinimizerWindowSize];
    MM2Settings::GapOpenDelete = options[OptionNames::GapOpenDelete];
    MM2Settings::GapOpenInsert = options[OptionNames::GapOpenInsert];
    MM2Settings::GapExtensionDelete = options[OptionNames::GapExtensionDelete];
    MM2Settings::GapExtensionInsert = options[OptionNames::GapExtensionInsert];
    MM2Settings::MatchScore = options[OptionNames::MatchScore];
    MM2Settings::MismatchPenalty = options[OptionNames::MismatchPenalty];
    MM2Settings::Zdrop = options[OptionNames::Zdrop];
    MM2Settings::ZdropInv = options[OptionNames::ZdropInv];
    MM2Settings::Bandwidth = options[OptionNames::Bandwidth];
    MM2Settings::MaxIntronLength = options[OptionNames::MaxIntronLength];
    MM2Settings::NonCanon = options[OptionNames::NonCanon];
    MM2Settings::NoSpliceFlank = options[OptionNames::NoSpliceFlank];

    int32_t requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    MM2Settings::NumThreads = ThreadCount(requestedNThreads);

    const std::map<std::string, AlignmentMode> alignModeMap{{"SUBREAD", AlignmentMode::SUBREADS},
                                                            {"ISOSEQ", AlignmentMode::ISOSEQ}};
    MM2Settings::AlignMode = alignModeMap.at(options[OptionNames::AlignModeOpt].get<std::string>());
}

int32_t AlignSettings::ThreadCount(int32_t n)
{
    const int32_t m = std::thread::hardware_concurrency();
    if (n <= 0) n = m + n;  // permit n <= 0 to subtract from max threads
    return std::max(1, std::min(m, n));
}

PacBio::CLI::Interface AlignSettings::CreateCLI()
{
    using Task = PacBio::CLI::ToolContract::Task;

    const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
    PacBio::CLI::Interface i{"pbmm2_align", "Align PacBio reads to reference sequences", version};

    // clang-format off
    i.AddGroup("Basic Options", {
        OptionNames::HelpOption,
        OptionNames::VersionOption,
        OptionNames::LogFile,
        OptionNames::LogLevelOption,
        OptionNames::NumThreads,
        OptionNames::ChunkSize,
    });

    i.AddGroup("Parameter Set Options", {
        OptionNames::AlignModeOpt,
    });

    i.AddGroup("General Parameter Override Options", {
        OptionNames::Kmer,
        OptionNames::MinimizerWindowSize,
        OptionNames::MatchScore,
        OptionNames::MismatchPenalty,
        OptionNames::GapOpenDelete,
        OptionNames::GapOpenInsert,
        OptionNames::GapExtensionDelete,
        OptionNames::GapExtensionInsert,
        OptionNames::Zdrop,
        OptionNames::ZdropInv,
        OptionNames::Bandwidth,
    });

    i.AddGroup("IsoSeq Parameter Override Options", {
        OptionNames::MaxIntronLength,
        OptionNames::NonCanon,
        OptionNames::NoSpliceFlank,
    });

    i.AddGroup("Read Group Options", {
        OptionNames::SampleName
    });

    i.AddGroup("Filter Options", {
        OptionNames::MinAccuracy,
        OptionNames::MinAlignmentLength,
        OptionNames::MedianFilter,
    });

    i.AddPositionalArguments({
        { "in.bam|xml", "Input BAM or DataSet XML", "<in.bam|xml>" },
        { "ref.fa|xml|mmi", "Reference FASTA, ReferenceSet XML, or Reference Index", "<ref.fa|xml|mmi>" },
        { "out.aligned.bam|xml", "Output BAM or DataSet XML", "[out.aligned.bam|xml]" }
    });

    const std::string id = "mapping.tasks.pbmm2_align";
    Task tcTask(id);
    tcTask.NumProcessors(Task::MAX_NPROC);

    tcTask.InputFileTypes({
        {
            "subread_set",
            "SubreadSet",
            "Subread DataSet or .bam file",
            "PacBio.DataSet.SubreadSet"
        },
        {
            "reference_set",
            "ReferenceSet",
            "ReferenceSet or .fasta file",
            "PacBio.DataSet.ReferenceSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "aligned_bam_output",
            "AlignmentSet",
            "AlignmentSet for output .bam file",
            "PacBio.DataSet.AlignmentSet",
            "pbmm2_output"
        }
    });

    CLI::ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}  // namespace minimap2
}  // namespace PacBio
