// Author: Armin Töpfer

#include <pbbam/BamHeader.h>
#include <pbbam/DataSet.h>
#include <pbbam/virtual/ZmwReadStitcher.h>
#include <pbcopper/logging/Logging.h>
#include <boost/algorithm/string.hpp>

#include <AlignSettings.h>
#include <InputOutputUX.h>
#include <Pbmm2Version.h>

#include <SampleNames.h>

namespace PacBio {
namespace minimap2 {
namespace {
static const std::string fallbackSampleName{"UnnamedSample"};
}
std::string SampleNames::SanitizeSampleName(const std::string& in)
{
    if (in.empty()) return fallbackSampleName;

    auto trimmed = boost::algorithm::trim_copy(in);
    if (trimmed.empty()) return fallbackSampleName;
    std::string sanitizedName;
    for (const char& c : trimmed) {
        if (c < '!' || c > '~')
            sanitizedName += '_';
        else
            sanitizedName += c;
    }
    return sanitizedName;
}

std::string SampleNames::SanitizeFileInfix(const std::string& in)
{
    std::string sanitizedName;
    for (const char& c : in) {
        if (c == '_' || c == '-' || c == '.' || (c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') ||
            (c >= 'a' && c <= 'z'))
            sanitizedName += c;
    }
    return sanitizedName;
}

MovieToSampleToInfix SampleNames::DetermineMovieToSampleToInfix(const BAM::DataSet& inFile)
{
    MovieToSampleToInfix movieNameToSampleAndInfix;
    const auto md = inFile.Metadata();
    if (md.HasChild("Collections")) {
        using DataSetElement = PacBio::BAM::internal::DataSetElement;

        DataSetElement collections = md.Child<DataSetElement>("Collections");
        for (const auto& collectionMetaData : collections.Children()) {
            if (!collectionMetaData->HasAttribute("Context")) {
                PBLOG_ERROR << "Cannot parse Context attribute of <CollectionMetadata> "
                               "element. Bailing on biosample parsing.";
                continue;
            }
            std::string movieName = collectionMetaData->Attribute("Context");
            std::string wellSampleName;
            std::string bioSampleName;
            const auto wellSample = collectionMetaData->Child<DataSetElement>("WellSample");
            if (wellSample.HasAttribute("Name")) wellSampleName = wellSample.Attribute("Name");
            if (wellSample.HasChild("BioSamples")) {
                const auto bioSamples = wellSample.Child<DataSetElement>("BioSamples");
                if (bioSamples.HasChild("BioSample")) {
                    const auto bioSample = bioSamples.Child<DataSetElement>("BioSample");
                    if (bioSample.HasAttribute("Name")) bioSampleName = bioSample.Attribute("Name");
                }
            }
            std::string finalName;
            if (!bioSampleName.empty())
                finalName = bioSampleName;
            else if (!wellSampleName.empty())
                finalName = wellSampleName;
            finalName = SanitizeSampleName(finalName);

            movieNameToSampleAndInfix[movieName] = {finalName, SanitizeFileInfix(finalName)};
        }
    }
    return movieNameToSampleAndInfix;
}

BAM::BamHeader SampleNames::GenerateBamHeader(const BAM::DataSet& inFile,
                                              const AlignSettings& settings, const UserIO& uio,
                                              const MovieToSampleToInfix& mtsti,
                                              std::string& fastxRgId)
{
    BAM::BamHeader hdr;
    if ((settings.HQRegion || settings.ZMW) && !uio.isAlignedInput) {
        BAM::ZmwReadStitcher reader(inFile);
        if (reader.HasNext()) {
            auto r = reader.Next();
            hdr = r.Header();
        }
    } else if (!uio.isFastaInput && !uio.isFastqInput) {
        const auto bamFiles = inFile.BamFiles();
        hdr = bamFiles.front().Header();
        for (size_t i = 1; i < bamFiles.size(); ++i)
            hdr += bamFiles.at(i).Header();
    } else {
        std::string rgString = settings.Rg;
        boost::replace_all(rgString, "\\t", "\t");
        if (rgString.empty()) rgString = "@RG\tID:default";
        BAM::ReadGroupInfo rg = BAM::ReadGroupInfo::FromSam(rgString);
        fastxRgId = rg.Id();
        if (rg.MovieName().empty()) rg.MovieName("default");
        hdr.AddReadGroup(rg);
    }
    if (uio.isAlignedInput) {
        hdr.ClearSequences();
    }

    bool performOverrideSampleName = !settings.SampleName.empty();
    std::string overridingSampleName = SampleNames::SanitizeSampleName(settings.SampleName);
    auto rgs = hdr.ReadGroups();
    hdr.ClearReadGroups();
    for (auto& rg : rgs) {
        if (performOverrideSampleName) {
            rg.Sample(overridingSampleName);
        } else if (rg.Sample().empty()) {
            if (uio.isFromXML || uio.isFromJson) {
                if (mtsti.find(rg.MovieName()) != mtsti.cend()) {
                    rg.Sample(mtsti.at(rg.MovieName()).first);
                } else {
                    PBLOG_INFO << "Cannot find biosample name for movie name " << rg.MovieName()
                               << "! Will use fallback.";
                    rg.Sample(SampleNames::SanitizeSampleName(""));
                }
            } else {
                rg.Sample(SampleNames::SanitizeSampleName(""));
            }
        }
        hdr.AddReadGroup(rg);
    }
    const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
    auto pg = BAM::ProgramInfo("pbmm2").Name("pbmm2").Version(version).CommandLine("pbmm2 " +
                                                                                   settings.CLI);
    if (!settings.TcOverrides.empty()) pg.CustomTags({{"or", settings.TcOverrides}});
    hdr.AddProgram(pg);
    return hdr;
}
}  // namespace minimap2
}  // namespace PacBio
