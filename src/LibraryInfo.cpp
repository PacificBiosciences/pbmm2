#include "LibraryGitHash.h"
#include "LibraryVersion.h"

#include <pbmm2/LibraryInfo.h>

#include <pbbam/LibraryInfo.h>

namespace PacBio {
namespace Pbmm2 {

Library::Info Minimap2LibraryInfo() { return {"minimap2", Pbmm2::Minimap2ReleaseVersion, ""}; }

Library::Bundle LibraryBundle()
{
    Library::Bundle bundle{LibraryInfo(), {}};
    bundle += Pbbam::LibraryBundle();
    bundle += Minimap2LibraryInfo();
    return bundle;
}

Library::Info LibraryInfo() { return {"pbmm2", Pbmm2::ReleaseVersion, Pbmm2::LibraryGitSha1}; }

}  // namespace Pbmm2
}  // namespace PacBio
