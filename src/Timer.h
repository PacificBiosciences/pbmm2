#pragma once

#include <chrono>
#include <string>

namespace PacBio {
class Timer
{
public:
    Timer();

    std::string ElapsedTime() const;
    void Restart();
    void Freeze();

public:
    static std::string ElapsedTimeFromSeconds(int64_t nanosecs);

private:
    std::chrono::time_point<std::chrono::steady_clock> tick_;
    bool frozen_ = false;
    std::chrono::time_point<std::chrono::steady_clock> tock_;
};
}  // namespace PacBio
