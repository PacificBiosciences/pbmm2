// Author: Armin Töpfer

#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/cli/HelpPrinter.h>
#include <pbcopper/cli/VersionPrinter.h>
#include <pbcopper/logging/Logging.h>

#include <Pbmm2Version.h>

#include "AlignSettings.h"
#include "AlignWorkflow.h"
#include "IndexSettings.h"
#include "IndexWorkflow.h"

struct Tool
{
    const std::string CodeName;
    const std::string HelpText;
    const std::string Example;
    const PacBio::CLI::Interface Interface;
    const PacBio::CLI::ResultsHandler Runner;
};

int main(int argc, char* argv[])
{
    std::map<std::string, Tool> tools;
    // clang-format off
    tools.insert(
        std::make_pair("align",
            Tool{"pbmm2_align",
                 "Align PacBio reads to reference sequences",
                 "ref.referenceset.xml movie.subreadset.xml ref.movie.alignmentset.xml",
                 PacBio::minimap2::AlignSettings::CreateCLI(),
                 &PacBio::minimap2::AlignWorkflow::Runner}));
    tools.insert(
        std::make_pair("index",
            Tool{"pbmm2_index",
                 "Index reference and store as .mmi file",
                 "ref.referenceset.xml ref.mmi",
                 PacBio::minimap2::IndexSettings::CreateCLI(),
                 &PacBio::minimap2::IndexWorkflow::Runner}));
    // clang-format on
    std::vector<std::string> toolOrder{"index", "align"};

    size_t toolNameWidth = 0;
    for (const auto& toolName_info : tools)
        toolNameWidth = std::max(toolNameWidth, toolName_info.first.size());

    if (argc == 2 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)) {
        std::cout << "pbmm2 " << PacBio::Pbmm2Version() << " (commit " << PacBio::Pbmm2GitSha1()
                  << ')' << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    const auto IsHelp = [](const char* c) {
        return strcmp(c, "-h") == 0 || strcmp(c, "--help") == 0;
    };
    const auto IsVersion = [](const char* c) { return strcmp(c, "--version") == 0; };

    if (argc == 1 || (argc == 2 && IsHelp(argv[1]))) {
        std::cout << "pbmm2 - minimap2 with native PacBio BAM support" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "  pbmm2 <tool>" << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -h, --help   Output this help." << std::endl;
        std::cout << "  --version    Output version info." << std::endl;
        std::cout << std::endl;
        std::cout << "Tools:" << std::endl;
        for (const auto& toolName : toolOrder) {
            std::cout << "    ";
            std::cout << std::setw(11) << std::left << toolName << tools.at(toolName).HelpText
                      << '\n';
        }
        std::cout << std::endl;
        std::cout << "Examples:" << std::endl;
        for (const auto& toolName_info : tools) {
            std::cout << "  pbmm2 ";
            std::cout << toolName_info.first << ' ' << toolName_info.second.Example << '\n';
        }
        std::cout << std::endl;
        std::cout << R"(Typical workflows:
  A. Generate index file for reference and reuse it to align reads
    $ pbmm2 index ref.fasta ref.mmi
    $ pbmm2 align ref.mmi movie.subreads.bam ref.movie.bam

  B. Align reads and sort on-the-fly, with 4 alignment and 2 sort threads
    $ pbmm2 align ref.fasta movie.subreads.bam ref.movie.bam --sort -j 4 -J 2

  C. Align reads, sort on-the-fly, and create PBI
    $ pbmm2 align ref.fasta movie.subreadset.xml ref.movie.alignmentset.xml --sort

  D. Omit output file and stream BAM output to stdout
    $ pbmm2 align hg38.mmi movie1.subreadset.xml | samtools sort > hg38.movie1.sorted.bam

  E. Align CCS fastq input and sort output
    $ pbmm2 align ref.fasta movie.Q20.fastq ref.movie.bam --sort --rg '@RG\tID:myid\tSM:mysample')"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    std::string toolName(argv[1]);

    const bool toolFound = tools.find(toolName) != tools.cend();
    if (!toolFound) {
        std::cerr << "ERROR: Unknown tool " << toolName << std::endl;
        return EXIT_FAILURE;
    }
    auto& selectedTool = tools.at(toolName);

    const auto PrintHelp = [&]() {
        std::stringstream ss;
        PacBio::CLI::HelpPrinter::Print(selectedTool.Interface, ss);
        std::string s = ss.str();
        const auto codeName = selectedTool.CodeName;
        size_t f = s.find(codeName);
        s.replace(f, codeName.length(), "pbmm2 " + toolName);
        std::cout << s;
    };
    const auto PrintVersion = [&]() {
        std::stringstream ss;
        PacBio::CLI::VersionPrinter::Print(selectedTool.Interface, ss);
        std::string s = ss.str();
        const auto codeName = selectedTool.CodeName;
        size_t f = s.find(codeName);
        s.replace(f, codeName.length(), "pbmm2 " + toolName);
        std::cout << s;
    };
    for (int i = 1; i < argc; ++i) {
        if (IsHelp(argv[i])) {
            PrintHelp();
            return EXIT_SUCCESS;
        }
        if (IsVersion(argv[i])) {
            PrintVersion();
            return EXIT_SUCCESS;
        }
    }
    if (argc == 2) {
        PrintHelp();
        return EXIT_SUCCESS;
    }

    try {
        return PacBio::CLI::Run(argc - 1, &argv[1], selectedTool.Interface, selectedTool.Runner);
    } catch (const std::runtime_error& e) {
        std::cerr << "ERROR: " << e.what();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
