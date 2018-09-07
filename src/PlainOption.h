// Author: Armin Töpfer

#pragma once

#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>

namespace PacBio {
namespace minimap2 {
class PlainOption
{
public:
    PlainOption(std::string id, std::vector<std::string> cliOptions, std::string name,
                std::string description, JSON::Json defaultValue,
                JSON::Json choices = JSON::Json(nullptr),
                CLI::OptionFlags flags = CLI::OptionFlags::DEFAULT)
        : id_(std::move(id))
        , cliOptions_(std::move(cliOptions))
        , name_(std::move(name))
        , description_(std::move(description))
        , defaultValue_(std::move(defaultValue))
        , choices_(std::move(choices))
        , flags_(std::move(flags))
    {}

    operator CLI::Option() const
    {
        return {id_, cliOptions_, description_, defaultValue_, choices_, flags_};
    }
    operator std::pair<std::string, std::string>() const { return std::make_pair(id_, name_); }
    operator std::string() const { return id_; }

    const JSON::Json& GetDefaultValue() const { return defaultValue_; }

private:
    std::string id_;
    std::vector<std::string> cliOptions_;
    std::string name_;
    std::string description_;
    JSON::Json defaultValue_;
    JSON::Json choices_ = JSON::Json(nullptr);
    CLI::OptionFlags flags_;
};
}  // namespace minimap2
}  // namespace PacBio
