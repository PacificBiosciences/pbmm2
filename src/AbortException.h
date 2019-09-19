// Author: Armin Töpfer

#pragma once

#include <stdexcept>

namespace PacBio {
namespace minimap2 {
class AbortException : public std::exception
{};
}  // namespace minimap2
}  // namespace PacBio
