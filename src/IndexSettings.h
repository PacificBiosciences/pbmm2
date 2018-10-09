// Author: Armin Töpfer

#pragma once

#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>
#include <pbmm2/MM2Settings.h>

#include "PlainOption.h"

namespace PacBio {
namespace minimap2 {
/// Contains user provided CLI configuration
struct IndexSettings : MM2Settings
{
    const std::string CLI;
    const std::vector<std::string> InputFiles;

    std::string LogFile;
    Logging::LogLevel LogLevel;

    int32_t BestN;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    IndexSettings(const PacBio::CLI::Results& options);

    int32_t ThreadCount(int32_t n);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI::Interface CreateCLI();
};
}  // namespace minimap2
}  // namespace PacBio
