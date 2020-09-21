// Author: Armin Töpfer
#include "SampleNames.h"

#include <fstream>

#include <pbbam/BamHeader.h>
#include <pbbam/DataSet.h>
#include <pbbam/virtual/ZmwReadStitcher.h>
#include <pbcopper/logging/Logging.h>
#include <boost/algorithm/string.hpp>

#include <pbmm2/Pbmm2Version.h>

#include "AbortException.h"
#include "AlignSettings.h"
#include "InputOutputUX.h"

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

MovieToSampleToInfix SampleNames::DetermineMovieToSampleToInfix(const UserIO& uio)
{
    MovieToSampleToInfix movieNameToSampleAndInfix;
    const auto FillForFile = [&](const std::string& f) {
        BAM::DataSet ds;
        try {
            ds = BAM::DataSet{f};
        } catch (...) {
            throw AbortException(UNKNOWN_FILE_TYPES);
        }

        // Check dataset's BAM header(s) for '@RG SM' tags.
        int namedSampleCount = 0;
        for (const auto& bamFile : ds.BamFiles()) {
            const auto& header = bamFile.Header();
            for (const auto& rg : header.ReadGroups()) {
                const auto movie = rg.MovieName();
                const auto sample = rg.Sample();
                if (!sample.empty()) {
                    movieNameToSampleAndInfix.emplace(
                        movie,
                        std::make_pair(SanitizeSampleName(sample), SanitizeFileInfix(sample)));
                    ++namedSampleCount;
                } else {
                    movieNameToSampleAndInfix.emplace(
                        movie, std::make_pair(SanitizeSampleName("UnnamedSample"),
                                              SanitizeFileInfix("UnnamedSample")));
                }
            }
        }

        const auto& md = ds.Metadata();
        const auto& biosamples = md.BioSamples();
        const auto biosampleCount = biosamples.Size();
        std::string nameFromMetadata;
        if (biosampleCount > 0) {
            if (biosampleCount > 1 && namedSampleCount == 0) {
                PBLOG_INFO << "Found more than 1 biosample, but read groups lack the SM tag - "
                              "using 'UnnamedSample'!";
                nameFromMetadata = "UnnamedSample";
            } else {
                nameFromMetadata = biosamples[0].Name();
            }
        }

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
                        if (bioSample.HasAttribute("Name"))
                            bioSampleName = bioSample.Attribute("Name");
                    }
                }
                std::string finalName;
                if (!nameFromMetadata.empty()) {
                    finalName = nameFromMetadata;
                } else if (!bioSampleName.empty())
                    finalName = bioSampleName;
                else if (!wellSampleName.empty())
                    finalName = wellSampleName;
                finalName = SanitizeSampleName(finalName);

                movieNameToSampleAndInfix[movieName] = {finalName, SanitizeFileInfix(finalName)};
            }
        }
    };
    if (uio.isFromJson && uio.isFromXML) {
        FillForFile(uio.unpackedFromJson);
    } else if (uio.isFromXML) {
        for (const auto& f : uio.inputFiles)
            FillForFile(f);
    }
    return movieNameToSampleAndInfix;
}

BAM::BamHeader SampleNames::GenerateBamHeader(const AlignSettings& settings, const UserIO& uio,
                                              const MovieToSampleToInfix& mtsti,
                                              std::string& fastxRgId)
{
    std::unique_ptr<BAM::BamHeader> hdr;
    if ((settings.HQRegion || settings.ZMW) && !uio.isAlignedInput) {
        if (uio.isFromFofn) {
            throw AbortException("Cannot combine --hqregion or --zmw with fofn input!");
        }
        BAM::ZmwReadStitcher reader(uio.inFile);
        if (reader.HasNext()) {
            auto r = reader.Next();
            hdr = std::make_unique<BAM::BamHeader>(r.Header().DeepCopy());
        }
    } else if (!uio.isFastaInput && !uio.isFastqInput) {
        const auto Fill = [&](const std::string& f) {
            BAM::DataSet ds{f};
            const auto bamFiles = ds.BamFiles();
            if (!hdr) {
                hdr = std::make_unique<BAM::BamHeader>(bamFiles.front().Header().DeepCopy());
                for (size_t i = 1; i < bamFiles.size(); ++i)
                    *hdr += bamFiles.at(i).Header();
            } else {
                for (size_t i = 0; i < bamFiles.size(); ++i)
                    *hdr += bamFiles.at(i).Header();
            }
        };
        if (uio.isFromJson) {
            Fill(uio.unpackedFromJson);
        } else {
            for (const auto& f : uio.inputFiles)
                Fill(f);
        };
    } else {
        hdr = std::make_unique<BAM::BamHeader>();
        std::string rgString = settings.Rg;
        boost::replace_all(rgString, "\\t", "\t");
        if (rgString.empty()) rgString = "@RG\tID:default";
        BAM::ReadGroupInfo rg = BAM::ReadGroupInfo::FromSam(rgString);
        fastxRgId = rg.Id();
        if (rg.MovieName().empty()) rg.MovieName("default");
        hdr->AddReadGroup(rg);
    }
    if (uio.isAlignedInput) {
        hdr->ClearSequences();
    }

    bool performOverrideSampleName = !settings.SampleName.empty();
    std::string overridingSampleName = SampleNames::SanitizeSampleName(settings.SampleName);
    auto rgs = hdr->ReadGroups();
    hdr->ClearReadGroups();
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
        hdr->AddReadGroup(rg);
    }
    const auto version = PacBio::Pbmm2Version() + " (commit " + PacBio::Pbmm2GitSha1() + ")";
    auto pg = BAM::ProgramInfo("pbmm2").Name("pbmm2").Version(version).CommandLine("pbmm2 " +
                                                                                   settings.CLI);
    hdr->AddProgram(pg);
    return hdr->DeepCopy();
}
}  // namespace minimap2
}  // namespace PacBio
