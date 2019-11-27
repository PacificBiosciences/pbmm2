// Author: Armin Töpfer

#pragma once

#include <stdexcept>
#include <string>
#include <utility>

namespace PacBio {
namespace minimap2 {

class AbortException : public std::exception
{
public:
    AbortException(std::string message) : message_{std::move(message)} {}

    const char* what() const noexcept override { return message_.c_str(); }

private:
    std::string message_;
};

}  // namespace minimap2
}  // namespace PacBio
