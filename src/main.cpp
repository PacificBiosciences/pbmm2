// Author: Armin Töpfer

#include <cstring>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/cli/HelpPrinter.h>
#include <pbcopper/cli/VersionPrinter.h>
#include <pbcopper/logging/Logging.h>

#include "Settings.h"
#include "Workflow.h"

int main(int argc, char* argv[])
{
    if (argc == 1) {
        PacBio::CLI::HelpPrinter::Print(PacBio::minimap2::Settings::CreateCLI(), std::cout);
        return EXIT_SUCCESS;
    }

    const auto IsHelp = [](const char* c) {
        return strcmp(c, "-h") == 0 || strcmp(c, "--help") == 0;
    };
    const auto IsVersion = [](const char* c) { return strcmp(c, "--version") == 0; };

    for (int i = 1; i < argc; ++i) {
        if (IsHelp(argv[i])) {
            PacBio::CLI::HelpPrinter::Print(PacBio::minimap2::Settings::CreateCLI(), std::cout);
            return EXIT_SUCCESS;
        }
        if (IsVersion(argv[i])) {
            PacBio::CLI::VersionPrinter::Print(PacBio::minimap2::Settings::CreateCLI(), std::cout);
            return EXIT_SUCCESS;
        }
    }

    try {
        return PacBio::CLI::Run(argc, argv, PacBio::minimap2::Settings::CreateCLI(),
                                &PacBio::minimap2::Workflow::Runner);
    } catch (const std::runtime_error& e) {
        PBLOG_FATAL << e.what();
    }
}