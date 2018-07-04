#include <iostream>
#include <stdexcept>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/cli/HelpPrinter.h>
#include <pbcopper/logging/Logging.h>

#include "Settings.h"
#include "Workflow.h"

int main(int argc, char* argv[])
{
    try {
        if (argc == 1)
            PacBio::CLI::HelpPrinter::Print(PacBio::minimap2::Settings::CreateCLI(), std::cout);
        else
            return PacBio::CLI::Run(argc, argv, PacBio::minimap2::Settings::CreateCLI(),
                                    &PacBio::minimap2::Workflow::Runner);
    } catch (const std::runtime_error& e) {
        PBLOG_FATAL << e.what();
    }
}