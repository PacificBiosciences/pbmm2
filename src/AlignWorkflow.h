#pragma once

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace minimap2 {
struct AlignWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};
}  // namespace minimap2
}  // namespace PacBio
