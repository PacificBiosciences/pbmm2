// Author: Armin Töpfer

#pragma once

#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>

#include "AlignmentMode.h"
#include "PlainOption.h"

namespace PacBio {
namespace minimap2 {
/// Contains user provided CLI configuration
struct AlignSettings
{
    const std::string CLI;
    const std::vector<std::string> InputFiles;
    int32_t NumThreads;

    // alignment output filters
    float MinAccuracy;
    int32_t MinAlignmentLength;

    bool Pbi;

    std::string LogFile;
    Logging::LogLevel LogLevel;

    const std::string SampleName;
    AlignmentMode AlignMode;
    int32_t BestN;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    AlignSettings(const PacBio::CLI::Results& options);

    int32_t ThreadCount(int32_t n);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI::Interface CreateCLI();
};
}  // namespace minimap2
}  // namespace PacBio
