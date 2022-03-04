#pragma once

#include <pbmm2/AlignmentMode.h>

#include <cstdint>
#include <string>

namespace PacBio {
namespace minimap2 {
struct MM2Settings
{
    AlignmentMode AlignMode = AlignmentMode::SUBREADS;
    int32_t NumThreads = -1;
    int32_t Kmer = -1;
    int32_t MinimizerWindowSize = -1;
    int32_t GapOpen1 = -1;
    int32_t GapOpen2 = -1;
    int32_t GapExtension1 = -1;
    int32_t GapExtension2 = -1;
    int32_t MatchScore = -1;
    int32_t MismatchPenalty = -1;
    int32_t Zdrop = -1;
    int32_t ZdropInv = -1;
    int32_t Bandwidth = -1;
    int32_t MaxIntronLength = -1;
    int32_t NonCanon = -1;
    int32_t MaxNumAlns = 0;
    int32_t MaxGap = -1;
    int32_t MaxSecondaryAlns = -1;
    bool NoSpliceFlank = false;
    bool DisableHPC = false;
    bool NoTrimming = false;
    bool ShortSACigar = false;
    float LongJoinFlankRatio = -1;
    std::string EnforcedMapping;
};
}  // namespace minimap2
}  // namespace PacBio
