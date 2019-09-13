// Author: Armin Töpfer

#pragma once

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace minimap2 {
struct IndexWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};
}  // namespace minimap2
}  // namespace PacBio
