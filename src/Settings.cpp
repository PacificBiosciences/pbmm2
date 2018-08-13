// Author: Armin Töpfer

#include <map>

#include <Pbmm2Version.h>

#include "Settings.h"

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
    "Minimum alignment accuracy.",
    CLI::Option::FloatType(0.75f)
};
const PlainOption MinAlignmentLength{
    "minalnlength",
    { "min-length" },
    "Minimum Alignment Length",
    "Minimum alignment length.",
    CLI::Option::IntType(50)
};
const PlainOption Pbi{
    "pbi",
    { "pbi" },
    "Generate PBI",
    "Generate a PBI file that is needed for SMRTLink, slower.",
    CLI::Option::BoolType()
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
    R"(Set alignment mode: "SUBREAD".)",
    CLI::Option::StringType("SUBREAD"),
    {"SUBREAD"}
};
const PlainOption BestN{
    "bestn",
    { "bestn" },
    "Max alignments",
    "Retain at most N alignments.",
    CLI::Option::IntType(5)
};
// clang-format on
}  // namespace OptionNames

Settings::Settings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , MinAccuracy(options[OptionNames::MinAccuracy])
    , MinAlignmentLength(options[OptionNames::MinAlignmentLength])
    , Pbi(options[OptionNames::Pbi])
    , LogFile{options[OptionNames::LogFile].get<decltype(LogFile)>()}
    , LogLevel{options.LogLevel()}
    , SampleName{options[OptionNames::SampleName].get<decltype(SampleName)>()}
    , BestN(options[OptionNames::BestN])
{
    int32_t requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    NumThreads = ThreadCount(requestedNThreads);

    const std::map<std::string, AlignmentMode> alignModeMap{{"SUBREAD", AlignmentMode::SUBREADS}};
    AlignMode = alignModeMap.at(options[OptionNames::AlignModeOpt].get<std::string>());
}

int32_t Settings::ThreadCount(int32_t n)
{
    const int32_t m = std::thread::hardware_concurrency();
    if (n <= 0) n = m + n;  // permit n <= 0 to subtract from max threads
    return std::max(1, std::min(m, n));
}

PacBio::CLI::Interface Settings::CreateCLI()
{
    using Task = PacBio::CLI::ToolContract::Task;

    const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
    PacBio::CLI::Interface i{"pbmm2", "minimap2 with native PacBio BAM support", version};

    // clang-format off
    i.AddGroup("Basic Options", {
        OptionNames::AlignModeOpt,
        OptionNames::HelpOption,
        OptionNames::VersionOption,
        OptionNames::LogFile,
        OptionNames::LogLevelOption,
        OptionNames::NumThreads
    });

    i.AddGroup("Read Group Options", {
        OptionNames::SampleName
    });

    i.AddGroup("Post-processing Options", {
        OptionNames::Pbi
    });

    i.AddGroup("Filter Options", {
        OptionNames::MinAccuracy,
        OptionNames::MinAlignmentLength,
        OptionNames::BestN
    });

    i.AddPositionalArguments({
        { "in.subreads.bam|xml", "Input BAM or DataSet XML", "<in.subreads.bam|xml>" },
        { "ref.fa|xml", "Reference FASTA or ReferenceSet XML", "<ref.fa|xml>" },
        { "out.aligned.bam|xml", "Output BAM or DataSet XML", "<out.aligned.bam|xml>" }
    });

    const std::string id = "mapping.tasks.pbmm2";
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