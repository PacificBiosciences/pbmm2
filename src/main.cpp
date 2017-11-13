#include <stdexcept>

#include <pbcopper/cli/CLI.h>

#include "Settings.h"
#include "Workflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI::Run(argc, argv, PacBio::minimap2::Settings::CreateCLI(),
                                &PacBio::minimap2::Workflow::Runner);
    } catch (const std::runtime_error& e) {
        std::cerr << "ERROR: " << e.what();
    }
}