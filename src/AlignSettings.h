// Author: Armin Töpfer

#pragma once

#include <string>
#include <vector>

#include <pbcopper/cli2/CLI.h>

#include "BamIndex.h"

#include <pbmm2/MM2Settings.h>

namespace PacBio {
namespace minimap2 {
/// Contains user provided CLI configuration
struct AlignSettings : MM2Settings
{
    const std::string CLI;
    const std::vector<std::string> InputFiles;

    // alignment output filters
    double MinPercConcordance;
    double MinPercIdentity;
    double MinPercIdentityGapComp;
    int32_t MinAlignmentLength;

    const std::string SampleName;
    int32_t ChunkSize;

    bool MedianFilter;

    bool Sort;
    int SortThreads;
    int64_t SortMemory;

    bool ZMW;
    bool HQRegion;

    bool Strip;
    bool SplitBySample;

    std::string Rg;

    bool CreatePbi;
    BamIndex BamIdx = BamIndex::NONE;
    bool OutputUnmapped;

    bool CompressSequenceHomopolymers;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    AlignSettings(const PacBio::CLI_v2::Results& options);

    int32_t ThreadCount(int32_t n);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI_v2::Interface CreateCLI();
};
}  // namespace minimap2
}  // namespace PacBio
