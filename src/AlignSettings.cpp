// Author: Armin Töpfer

#include <unistd.h>
#include <map>

#include <boost/algorithm/string.hpp>

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
    { "j", "alignment-threads" },
    "Number of Threads",
    "Number of threads used for alignment, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption MinPercConcordance{
    "min_perc_concordance",
    { "c", "min-concordance-perc" },
    "Minimum Concordance (%)",
    "Minimum alignment concordance in percent.",
    CLI::Option::IntType(70)
};
const PlainOption MinAlignmentLength{
    "minalnlength",
    { "l", "min-length" },
    "Minimum Length (bp)",
    "Minimum mapped read length in basepair.",
    CLI::Option::IntType(50)
};
const PlainOption SampleName{
    "biosample_name",
    { "sample" },
    "Sample Name",
    "Sample name for all read groups. Defaults, in order of precedence: SM field in input read group, biosample name, well sample name, \"UnnamedSample\".",
    CLI::Option::StringType()
};
const PlainOption AlignModeOpt{
    "align_mode",
    { "preset" },
    "Alignment mode",
    "Set alignment mode:\n  - \"SUBREAD\" -k 19 -w 10 -o 5 -O 56 -e 4 -E 1 -A 2 -B 5 -z 400 -Z 50 -r 2000 -L 0.5\n"
    "  - \"CCS\" -k 19 -w 10 -u -o 5 -O 56 -e 4 -E 1 -A 2 -B 5 -z 400 -Z 50 -r 2000 -L 0.5\n"
    "  - \"ISOSEQ\" -k 15 -w 5 -u -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -C 5 -r 200000 -G 200000 -L 0.5\n"
    "  - \"UNROLLED\" -k 15 -w 15 -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 2000 -L 0.5\n"
    "Default",
    CLI::Option::StringType("SUBREAD"),
    {"SUBREAD", "CCS", "ISOSEQ", "UNROLLED"}
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
const PlainOption GapOpen1{
    "gap_open_1",
    { "o", "gap-open-1" },
    "Gap Open Penalty 1",
    "Gap open penalty 1.",
    CLI::Option::IntType(-1)
};
const PlainOption GapOpen2{
    "gap_open_2",
    { "O", "gap-open-2" },
    "Gap Open Penalty 2",
    "Gap open penalty 2.",
    CLI::Option::IntType(-1)
};
const PlainOption GapExtension1{
    "gap_extension_1",
    { "e", "gap-extend-1" },
    "Gap Extension Penalty 1",
    "Gap extension penalty 1.",
    CLI::Option::IntType(-1)
};
const PlainOption GapExtension2{
    "gap_extension_ins",
    { "E", "gap-extend-2" },
    "Gap Extension Penalty 2",
    "Gap extension penalty 2.",
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
    "Max Intron Length (effective in ISOSEQ preset; changing bandwidth)",
    "Max intron length (changes -r).",
    CLI::Option::IntType(-1)
};
const PlainOption NonCanon{
    "non_canon",
    { "C" },
    "Cost For a Non-Canonical GT-AG Splicing (effective in ISOSEQ preset)",
    "Cost for a non-canonical GT-AG splicing.",
    CLI::Option::IntType(-1)
};
const PlainOption NoSpliceFlank{
    "no_splice_flank",
    { "no-splice-flank" },
    "Do Not Prefer Splice Flanks GT-AG (effective in ISOSEQ preset)",
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
const PlainOption Sort{
    "sort",
    { "sort" },
    "Generate sorted BAM file",
    "Generate sorted BAM file.",
    CLI::Option::BoolType(false)
};
const PlainOption Pbi{
    "pbi",
    { "pbi" },
    "Generate PBI file",
    "Generate PBI file, only works with --sort.",
    CLI::Option::BoolType(false)
};
const PlainOption SortThreadsTC{
    "sort_threads_tc",
    { "sort-threads-perc" },
    "Percentage of threads used for sorting",
    "Percentage of threads used exclusively for sorting (absolute number of sort threads is capped at 8).",
    CLI::Option::IntType(25),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};
const PlainOption SortThreads{
    "sort_threads",
    { "J", "sort-threads" },
    "Number of threads used for sorting",
    "Number of threads used for sorting; 0 means 25% of -j, maximum 8.",
    CLI::Option::IntType(0)
};
const PlainOption SortMemory{
    "sort_memory",
    { "m", "sort-memory" },
    "Memory per thread for sorting",
    "Memory per thread for sorting.",
    CLI::Option::StringType("768M")
};
const PlainOption SortMemoryTC{
    "sort_memory_tc",
    { "sort-memory-tc" },
    "Memory per thread for sorting",
    "Memory per thread for sorting.",
    CLI::Option::StringType("4G"),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};
const PlainOption DisableHPC{
    "disable_hpc",
    { "u", "no-kmer-compression" },
    "Disable Homopolymer-Compressed seeding",
    "Disable homopolymer-compressed k-mer (compression is activate for SUBREAD & UNROLLED presets).",
    CLI::Option::BoolType(false)
};
const PlainOption ZMW{
    "zmw_mode",
    { "zmw" },
    "Process ZMW Reads",
    "Process ZMW Reads, subreadset.xml input required (activates UNROLLED preset).",
    CLI::Option::BoolType(false)
};
const PlainOption HQRegion{
    "hq_mode",
    { "hqregion" },
    "Process HQ Regions",
    "Process HQ region of each ZMW, subreadset.xml input required (activates UNROLLED preset).",
    CLI::Option::BoolType(false)
};
const PlainOption Strip{
    "strip",
    { "strip" },
    "Strip Base Tags",
    "Remove all kinetic and extra QV tags. Output cannot be polished.",
    CLI::Option::BoolType(false)
};
const PlainOption SplitBySample{
    "split_by_sample",
    { "split-by-sample" },
    "Split by Sample",
    "One output BAM per sample.",
    CLI::Option::BoolType(false)
};
const PlainOption Rg{
    "rg",
    { "rg" },
    "Read group",
    "Read group header line such as '@RG\\tID:xyz\\tSM:abc'. Only for FASTA/Q inputs.",
    CLI::Option::StringType()
};
const PlainOption CreatePbi{
    "create_pbi",
    { "pbi" },
    "Generate PBI for BAM only output",
    "Generate PBI for BAM only output.",
    CLI::Option::BoolType(false),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};
const PlainOption LongJoinFlankRatio{
    "long_join_flank_ratio",
    { "L", "lj-min-ratio" },
    "Long Join Flank Ratio",
    "Long join flank ratio.",
    CLI::Option::FloatType(-1)
};
const PlainOption NoBAI{
    "no_bai",
    { "no-bai" },
    "Omit BAI Generation",
    "Omit BAI generation for sorted output.",
    CLI::Option::BoolType(false)
};
const PlainOption NoTrimming{
    "disable_repeated_matches_trimming",
    { "no-rmt" },
    "No Repeated Matches Trimming",
    "Disable repeated matches trimming.",
    CLI::Option::BoolType(false),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};
// clang-format on
}  // namespace OptionNames

AlignSettings::AlignSettings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , IsFromRTC(options.IsFromRTC())
    , MinPercConcordance(options[OptionNames::MinPercConcordance])
    , MinAlignmentLength(options[OptionNames::MinAlignmentLength])
    , LogFile{options[OptionNames::LogFile].get<decltype(LogFile)>()}
    , LogLevel{options.LogLevel()}
    , SampleName{options[OptionNames::SampleName].get<decltype(SampleName)>()}
    , ChunkSize(options[OptionNames::ChunkSize])
    , MedianFilter(options[OptionNames::MedianFilter])
    , Sort(options[OptionNames::Sort])
    , ZMW(options[OptionNames::ZMW])
    , HQRegion(options[OptionNames::HQRegion])
    , Strip(options[OptionNames::Strip])
    , SplitBySample(options[OptionNames::SplitBySample])
    , Rg(options[OptionNames::Rg].get<decltype(Rg)>())
    , CreatePbi(options[OptionNames::CreatePbi])
    , NoBAI(options[OptionNames::NoBAI])
{
    MM2Settings::Kmer = options[OptionNames::Kmer];
    MM2Settings::MinimizerWindowSize = options[OptionNames::MinimizerWindowSize];
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
    MM2Settings::DisableHPC = options[OptionNames::DisableHPC];
    MM2Settings::LongJoinFlankRatio = options[OptionNames::LongJoinFlankRatio];
    MM2Settings::NoTrimming = options[OptionNames::NoTrimming];

    int numAvailableCores = std::thread::hardware_concurrency();
    int32_t requestedNThreads;
    std::string requestedMemory;
    int sortThreadPerc = 25;

    if (IsFromRTC) {
        Sort = true;
        requestedNThreads = options.NumProcessors();
        requestedMemory = options[OptionNames::SortMemoryTC].get<std::string>();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
        requestedMemory = options[OptionNames::SortMemory].get<std::string>();
    }
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
    SortMemory = PlainOption::SizeStringToInt(requestedMemory);

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
        int64_t maxMem = SortMemory * SortThreads;
        float maxMemSortFloat;
        std::string maxMemSortSuffix;
        MemoryToHumanReadable(SortMemory * SortThreads, &maxMemSortFloat, &maxMemSortSuffix);
        PBLOG_INFO << "Using " << MM2Settings::NumThreads << " threads for alignments, "
                   << SortThreads << " threads for sorting, and " << maxMemSortFloat
                   << maxMemSortSuffix << " bytes RAM for sorting.";

        auto pages = sysconf(_SC_PHYS_PAGES);
        auto page_size = sysconf(_SC_PAGE_SIZE);
        auto availableMemory = pages * page_size;

        float availFloat;
        std::string availSuffix;
        MemoryToHumanReadable(availableMemory, &availFloat, &availSuffix);

        if (maxMem > availableMemory) {
            PBLOG_FATAL << "Trying to allocate more memory for sorting (" << maxMemSortFloat
                        << maxMemSortSuffix << ") than system-wide available (" << availFloat
                        << availSuffix << ")";
            std::exit(EXIT_FAILURE);
        }
    } else {
        PBLOG_INFO << "Using " << MM2Settings::NumThreads << " threads for alignments.";
    }

    const std::map<std::string, AlignmentMode> alignModeMap{{"SUBREAD", AlignmentMode::SUBREADS},
                                                            {"ISOSEQ", AlignmentMode::ISOSEQ},
                                                            {"CCS", AlignmentMode::CCS},
                                                            {"UNROLLED", AlignmentMode::UNROLLED}};

    MM2Settings::AlignMode = alignModeMap.at(options[OptionNames::AlignModeOpt].get<std::string>());
    int inputFilterCounts = ZMW + MedianFilter + HQRegion;
    if (inputFilterCounts > 1) {
        PBLOG_FATAL << "Options --zmw, --hqregion and --median-filter are mutually exclusive.";
        std::exit(EXIT_FAILURE);
    }
    if (ZMW || HQRegion) {
        if (ChunkSize != 100)
            PBLOG_WARN << "Cannot change --chunk-size in --zmw/--hqregion mode. Parameters "
                          "--chunk-size is forced to 1.";
        ChunkSize = 1;
        MM2Settings::AlignMode = AlignmentMode::UNROLLED;
    }

    if (!Rg.empty() && !boost::contains(Rg, "ID") && !boost::starts_with(Rg, "@RG\t")) {
        PBLOG_FATAL << "Invalid @RG line. Missing ID field. Please provide following "
                       "format: '@RG\\tID:xyz\\tSM:abc'";
        std::exit(EXIT_FAILURE);
    }
    if (MM2Settings::LongJoinFlankRatio > 1) {
        PBLOG_FATAL << "Option -L,--lj-min-ratio has to be between a ratio betweem 0 and 1.";
        std::exit(EXIT_FAILURE);
    }

    if (!Sort && NoBAI) {
        PBLOG_WARN << "Option --no-bai has no effect without option --sort!";
    }
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
        OptionNames::ChunkSize,
        OptionNames::NoTrimming,

        // hidden
        OptionNames::SortMemoryTC,
        OptionNames::CreatePbi,
    });

    i.AddGroup("Sorting Options", {
        OptionNames::Sort,
        OptionNames::SortMemory,
    });

    i.AddGroup("Threading Options", {
        OptionNames::NumThreads,
        OptionNames::SortThreads,
        OptionNames::SortThreadsTC,
    });

    i.AddGroup("Parameter Set Options", {
        OptionNames::AlignModeOpt,
    });

    i.AddGroup("General Parameter Override Options", {
        OptionNames::Kmer,
        OptionNames::MinimizerWindowSize,
        OptionNames::DisableHPC,
        OptionNames::MatchScore,
        OptionNames::MismatchPenalty,
        OptionNames::Zdrop,
        OptionNames::ZdropInv,
        OptionNames::Bandwidth,
    });

    i.AddGroup("Gap Parameter Override Options (a k-long gap costs min{o+k*e,O+k*E})", {
        OptionNames::GapOpen1,
        OptionNames::GapOpen2,
        OptionNames::GapExtension1,
        OptionNames::GapExtension2,
        OptionNames::LongJoinFlankRatio,
    });

    i.AddGroup("IsoSeq Parameter Override Options", {
        OptionNames::MaxIntronLength,
        OptionNames::NonCanon,
        OptionNames::NoSpliceFlank,
    });

    i.AddGroup("Read Group Options", {
        OptionNames::SampleName,
        OptionNames::Rg,
    });

    i.AddGroup("Output Options", {
        OptionNames::MinPercConcordance,
        OptionNames::MinAlignmentLength,
        OptionNames::Strip,
        OptionNames::SplitBySample,
        OptionNames::NoBAI,
    });

    i.AddGroup("Input Manipulation Options (mutually exclusive)", {
        OptionNames::MedianFilter,
        OptionNames::ZMW,
        OptionNames::HQRegion,
    });

    i.AddGroup("Output File Options", {

    });

    i.AddPositionalArguments({
        { "ref.fa|xml|mmi", "Reference FASTA, ReferenceSet XML, or Reference Index", "<ref.fa|xml|mmi>" },
        { "in.bam|xml|fa|fq", "Input BAM, DataSet XML, FASTA, or FASTQ", "<in.bam|xml|fa|fq>" },
        { "out.aligned.bam|xml", "Output BAM or DataSet XML", "[out.aligned.bam|xml]" }
    });

    const std::string id = "mapping.tasks.pbmm2_align";
    Task tcTask(id);
    tcTask.NumProcessors(Task::MAX_NPROC);
    tcTask.AddOption(OptionNames::MinPercConcordance);
    tcTask.AddOption(OptionNames::MinAlignmentLength);
    tcTask.AddOption(OptionNames::SortMemoryTC);
    tcTask.AddOption(OptionNames::SampleName);
    tcTask.AddOption(OptionNames::ZMW);
    tcTask.AddOption(OptionNames::HQRegion);
    tcTask.AddOption(OptionNames::MedianFilter);
    tcTask.AddOption(OptionNames::Strip);
    tcTask.AddOption(OptionNames::SplitBySample);

    tcTask.InputFileTypes({
        {
            "datastore_input",
            "Datastore",
            "Datastore containing ONE dataset",
            "PacBio.FileTypes.json"
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
            "datastore_output",
            "Datastore",
            "Datastore containing one dataset",
            "PacBio.FileTypes.json",
            "out"
        },
    });

    CLI::ToolContract::Driver driver;
    driver.Exe("pbmm2 align --resolved-tool-contract");
    CLI::ToolContract::Config tcConfig(tcTask, driver);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}  // namespace minimap2
}  // namespace PacBio
