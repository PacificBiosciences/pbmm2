// Author: Armin Töpfer

#pragma once

#include <cstdint>

#include "AlignmentMode.h"

namespace PacBio {
namespace minimap2 {
struct MM2Settings
{
    AlignmentMode AlignMode;
    int32_t NumThreads = -1;
    int32_t Kmer = -1;
    int32_t MinimizerWindowSize = -1;
    int32_t GapOpenDelete = -1;
    int32_t GapOpenInsert = -1;
    int32_t GapExtensionDelete = -1;
    int32_t GapExtensionInsert = -1;
    int32_t MatchScore = -1;
    int32_t MismatchPenalty = -1;
    int32_t Zdrop = -1;
    int32_t ZdropInv = -1;
    int32_t Bandwidth = -1;
    int32_t MaxIntronLength = -1;
    int32_t NonCanon = -1;
    bool NoSpliceFlank = false;
};
}  // namespace minimap2
}  // namespace PacBio
