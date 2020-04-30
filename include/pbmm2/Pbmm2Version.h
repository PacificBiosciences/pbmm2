// Author: Armin Töpfer

#pragma once

#include <string>

namespace PacBio {

std::string Pbmm2GitSha1();
std::string Pbmm2Version();

inline std::string Pbmm2FormattedVersion()
{
    return Pbmm2Version() + " (commit " + Pbmm2GitSha1() + ")";
}

std::string Minimap2Version();

}  // namespace PacBio
