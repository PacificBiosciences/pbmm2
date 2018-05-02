#include "Settings.h"

namespace PacBio {
namespace minimap2 {
namespace OptionNames {
// clang-format off
const PlainOption NumThreads{
    "numthreads",
    { "j", "num-threads" },
    "Number of Threads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption Kmer{
    "kmer",
    { "k", "kmer" },
    "K-mer Length",
    "The K in K-mer.",
    CLI::Option::IntType(19)
};
const PlainOption Window{
    "window",
    { "w", "window" },
    "Minimizer Window Size",
    "Minimizer window size.",
    CLI::Option::IntType(10)
};
const PlainOption NoHPC{
    "nohpc",
    { "no-hpc" },
    "Disable Homopolymer Compression",
    "Disable homopolymer compression.",
    CLI::Option::BoolType(false)
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
const PlainOption NoPbi{
    "nopbi",
    { "no-pbi" },
    "No PBI",
    "Do not generate a PBI file that is needed for SMRTLink.",
    CLI::Option::BoolType()
};
// clang-format on
}  // namespace OptionNames

Settings::Settings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , Kmer(options[OptionNames::Kmer])
    , Window(options[OptionNames::Window])
    , NoHPC(options[OptionNames::NoHPC])
    , MinAccuracy(options[OptionNames::MinAccuracy])
    , MinAlignmentLength(options[OptionNames::MinAlignmentLength])
    , NoPbi(true)
{
    int requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    NumThreads = ThreadCount(requestedNThreads);
}

int Settings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();
    if (n <= 0) n = m + n;  // permit n <= 0 to subtract from max threads
    return std::max(1, std::min(m, n));
}

PacBio::CLI::Interface Settings::CreateCLI()
{
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"pbmm2", "minimap2 with native PacBio BAM support", "0.3.0"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddOptions({
        OptionNames::NumThreads, OptionNames::Kmer, OptionNames::Window, OptionNames::NoHPC,
        OptionNames::MinAccuracy, OptionNames::MinAlignmentLength
        // OptionNames::NoPbi
    });

    // clang-format off
    i.AddPositionalArguments({
        { "INPUT", "Source BAM or DataSet XML", "INPUT" },
        { "REFERENCE", "Reference FASTA or ReferenceSet XML", "REFERENCE" },
        { "OUTPUT", "Output BAM or DataSet XML", "OUTPUT" }
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