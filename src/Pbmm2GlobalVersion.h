#pragma once

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/MultiToolInterface.h>

namespace PacBio {
namespace Pbmm2 {
void PrintPbmm2VersionMulti(const CLI_v2::MultiToolInterface& interface);
void PrintPbmm2VersionSingle(const CLI_v2::Interface& interface);
}  // namespace Pbmm2
}  // namespace PacBio
