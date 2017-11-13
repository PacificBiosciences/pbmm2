#pragma once

#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>

namespace PacBio {
namespace minimap2 {
struct PlainOption
{
    std::string id;
    std::vector<std::string> cliOptions;
    std::string name;
    std::string description;
    JSON::Json defaultValue;
    JSON::Json choices = JSON::Json(nullptr);
    CLI::OptionFlags flags;

    PlainOption(const std::string id, const std::vector<std::string> cliOptions,
                const std::string name, const std::string description,
                const JSON::Json defaultValue, const JSON::Json choices = JSON::Json(nullptr),
                const CLI::OptionFlags flags = CLI::OptionFlags::DEFAULT)
        : id(std::move(id))
        , cliOptions(std::move(cliOptions))
        , name(std::move(name))
        , description(std::move(description))
        , defaultValue(std::move(defaultValue))
        , choices(std::move(choices))
        , flags(std::move(flags))
    {
    }

    operator CLI::Option() const
    {
        return {id, cliOptions, description, defaultValue, choices, flags};
    }
    operator std::pair<std::string, std::string>() const { return std::make_pair(id, name); }
    operator std::string() const { return id; }
};

/// Contains user provided CLI configuration for lima_raw
struct Settings
{
    const std::string CLI;
    const std::vector<std::string> InputFiles;
    int NumThreads;
    bool NoPbi;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    Settings(const PacBio::CLI::Results& options);

    size_t ThreadCount(int n);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI::Interface CreateCLI();
};
}
}  // ::PacBio::minimap2