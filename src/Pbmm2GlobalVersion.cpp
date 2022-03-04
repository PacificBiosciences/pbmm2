#include "Pbmm2GlobalVersion.h"

#include <pbmm2/LibraryInfo.h>

#include <htslib/hts.h>
#include <pbbam/PbbamVersion.h>
#include <pbcopper/utility/PbcopperVersion.h>
#include <zlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/version.hpp>

#include <iostream>

namespace PacBio {
namespace Pbmm2 {
namespace {
void PrintPbmm2Version(const std::string& applicationName, const std::string& applicationVersion)
{
    const std::string pbmm2Version = []() {
        return Pbmm2::LibraryInfo().Release + " (commit " + Pbmm2::LibraryInfo().GitSha1 + ')';
    }();
    const std::string mm2Version = []() { return Pbmm2::Minimap2LibraryInfo().Release; }();
    const std::string pbbamVersion = []() { return BAM::LibraryFormattedVersion(); }();
    const std::string pbcopperVersion = []() {
        return Utility::LibraryVersionString() + " (commit " + Utility::LibraryGitSha1String() +
               ')';
    }();
    const std::string boostVersion = []() {
        std::string v = BOOST_LIB_VERSION;
        boost::replace_all(v, "_", ".");
        return v;
    }();
    const std::string htslibVersion = []() { return std::string{hts_version()}; }();
    const std::string zlibVersion = []() { return std::string{ZLIB_VERSION}; }();

    std::cout << applicationName << " " << applicationVersion << '\n';
    std::cout << '\n';
    std::cout << "Using:\n";
    std::cout << "  pbmm2    : " << pbmm2Version << '\n';
    std::cout << "  pbbam    : " << pbbamVersion << '\n';
    std::cout << "  pbcopper : " << pbcopperVersion << '\n';
    std::cout << "  boost    : " << boostVersion << '\n';
    std::cout << "  htslib   : " << htslibVersion << '\n';
    std::cout << "  minimap2 : " << mm2Version << '\n';
    std::cout << "  zlib     : " << zlibVersion << '\n';
}
}  // namespace
void PrintPbmm2VersionMulti(const CLI_v2::MultiToolInterface& interface)
{
    PrintPbmm2Version(interface.ApplicationName(), interface.ApplicationVersion());
}
void PrintPbmm2VersionSingle(const CLI_v2::Interface& interface)
{
    PrintPbmm2Version(interface.ApplicationName(), interface.ApplicationVersion());
}
}  // namespace Pbmm2
}  // namespace PacBio
