// Author: Armin Töpfer

#include <map>

#include <Pbmm2Version.h>

#include "IndexSettings.h"

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
const PlainOption AlignModeOpt{
    "align_mode",
    { "preset" },
    "Alignment mode",
    "Set alignment mode:\n  - \"SUBREAD\" -k 19 -w 10\n  - \"CCS\"  -k 19 -w 10 --no-hpc\n  - \"ISOSEQ\"  -k 15 -w 5 --no-hpc\nDefault",
    CLI::Option::StringType("SUBREAD"),
    {"SUBREAD", "CCS", "ISOSEQ"}
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
const PlainOption DisableHPC{
    "disable_hpc",
    { "no-hpc" },
    "Disable Homopolymer-Compressed seeding",
    "Disable homopolymer-compressed k-mer (hpc is only activated for SUBREAD mode).",
    CLI::Option::BoolType(false)
};
// clang-format on
}  // namespace OptionNames

IndexSettings::IndexSettings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , LogFile{options[OptionNames::LogFile].get<decltype(LogFile)>()}
    , LogLevel{options.LogLevel()}
{
    MM2Settings::Kmer = options[OptionNames::Kmer];
    MM2Settings::MinimizerWindowSize = options[OptionNames::MinimizerWindowSize];
    MM2Settings::DisableHPC = options[OptionNames::DisableHPC];

    int32_t requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    MM2Settings::NumThreads = ThreadCount(requestedNThreads);

    const std::map<std::string, AlignmentMode> alignModeMap{{"SUBREAD", AlignmentMode::SUBREADS},
                                                            {"CCS", AlignmentMode::CCS},
                                                            {"ISOSEQ", AlignmentMode::ISOSEQ}};
    MM2Settings::AlignMode = alignModeMap.at(options[OptionNames::AlignModeOpt].get<std::string>());
}

int32_t IndexSettings::ThreadCount(int32_t n)
{
    const int32_t m = std::thread::hardware_concurrency();
    if (n <= 0) n = m + n;  // permit n <= 0 to subtract from max threads
    return std::max(1, std::min(m, n));
}

PacBio::CLI::Interface IndexSettings::CreateCLI()
{
    using Task = PacBio::CLI::ToolContract::Task;

    const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
    PacBio::CLI::Interface i{"pbmm2_index", "Index reference and store as .mmi file", version};

    // clang-format off
    i.AddGroup("Basic Options", {
        OptionNames::HelpOption,
        OptionNames::VersionOption,
        OptionNames::LogFile,
        OptionNames::LogLevelOption,
        OptionNames::NumThreads,
    });

    i.AddGroup("Parameter Set Option", {
        OptionNames::AlignModeOpt,
    });

    i.AddGroup("Parameter Override Options", {
        OptionNames::Kmer,
        OptionNames::MinimizerWindowSize,
        OptionNames::DisableHPC
    });

    i.AddPositionalArguments({
        { "ref.fa|xml", "Reference FASTA, ReferenceSet XML", "<ref.fa|xml>" },
        { "out.mmi", "Output Reference Index", "<out.mmi>" }
    });

    const std::string id = "mapping.tasks.pbmm2_index";
    Task tcTask(id);
    tcTask.NumProcessors(Task::MAX_NPROC);

    tcTask.InputFileTypes({
        {
            "reference_set",
            "ReferenceSet",
            "ReferenceSet or .fasta file",
            "PacBio.DataSet.ReferenceSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "reference_index_output",
            "Minimap2IndexSet",
            "Minimap2IndexSet for output .mmi file",
            "PacBio.DataSet.Minimap2IndexSet",
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
