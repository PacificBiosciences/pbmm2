#pragma once

#include <pbmm2/AlignmentMode.h>
#include "AlignSettings.h"
#include "InputOutputUX.h"

#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>

#include <map>
#include <string>
#include <vector>

namespace PacBio {
namespace minimap2 {

using MovieToSampleToInfix = std::map<std::string, std::pair<std::string, std::string>>;

class SampleNames
{
public:
    static std::string SanitizeSampleName(const std::string& in);
    static std::string SanitizeFileInfix(const std::string& in);
    static MovieToSampleToInfix DetermineMovieToSampleToInfix(const UserIO& uio);
    static BAM::BamHeader GenerateBamHeader(const AlignSettings& settings, const UserIO& uio,
                                            const MovieToSampleToInfix& mtsti,
                                            std::string& fastxRgId);
};
}  // namespace minimap2
}  // namespace PacBio
