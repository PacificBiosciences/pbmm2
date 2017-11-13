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

size_t Settings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();
    if (n <= 0) n = m + n;  // permit n <= 0 to subtract from max threads
    return std::max(1, std::min(m, n));
}

PacBio::CLI::Interface Settings::CreateCLI()
{
    using Option = PacBio::CLI::Option;
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"pbmm2", "minimap2 with native PacBio BAM support",
                             "0.0.1"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddOptions({
        OptionNames::NumThreads
        // OptionNames::NoPbi
    });

    // clang-format off
    i.AddPositionalArguments({
        { "input", "Source BAM or DATASET", "INPUT" },
        { "reference", "FASTA", "REFERENCE" },
        { "output", "Output BAM", "OUTPUT" }
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
} // namespace minimap2
} // namespace PacBio