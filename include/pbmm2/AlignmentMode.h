#pragma once

namespace PacBio {
namespace minimap2 {
enum class AlignmentMode
{
    SUBREADS = 0,
    ISOSEQ,
    CCS,
    UNROLLED
};
}
}  // namespace PacBio
