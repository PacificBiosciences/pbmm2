#pragma once

#include <pbcopper/library/Bundle.h>
#include <pbcopper/library/Info.h>

namespace PacBio {
namespace Pbmm2 {

///
/// \return pbmm2 library info (e.g. name, version)
///
Library::Info LibraryInfo();

///
/// \returns bundle of pbmm2 library info, plus its dependencies
///
Library::Bundle LibraryBundle();

///
/// \return minimap2 library info (pbmm2 dependency)
///
Library::Info Minimap2LibraryInfo();

}  // namespace Pbmm2
}  // namespace PacBio
