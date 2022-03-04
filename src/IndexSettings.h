#pragma once

#include <pbmm2/MM2Settings.h>

#include <pbcopper/cli2/CLI.h>

#include <string>
#include <vector>

namespace PacBio {
namespace minimap2 {
/// Contains user provided CLI configuration
struct IndexSettings : MM2Settings
{
    const std::string CLI;
    const std::vector<std::string> InputFiles;
    int32_t BestN;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    IndexSettings(const PacBio::CLI_v2::Results& options);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI_v2::Interface CreateCLI();
};
}  // namespace minimap2
}  // namespace PacBio
