
# Get static libraries
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Boost
if(NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

# pbcopper
if (NOT pbcopper_INCLUDE_DIRS OR
    NOT pbcopper_LIBRARIES)
    set(pbcopper_build_tests OFF CACHE INTERNAL "" FORCE)
    set(pbcopper_build_docs OFF CACHE INTERNAL "" FORCE)
    set(pbcopper_build_examples OFF CACHE INTERNAL "" FORCE)
    add_subdirectory(${PBMM2_ThirdPartyDir}/pbcopper external/pbcopper/build)
endif()

# Threads
if (NOT Threads)
    find_package(Threads REQUIRED)
endif()

# pbbam
if (NOT PacBioBAM_INCLUDE_DIRS OR
    NOT PacBioBAM_LIBRARIES)
    set(PacBioBAM_build_docs    OFF CACHE INTERNAL "" FORCE)
    set(PacBioBAM_build_tests   OFF CACHE INTERNAL "" FORCE)
    set(PacBioBAM_build_tools   OFF CACHE INTERNAL "" FORCE)
    add_subdirectory(${PBMM2_ThirdPartyDir}/pbbam external/pbbam/build)
endif()

# minimap2
if (NOT minimap2_INCLUDE_DIRS OR
    NOT minimap2_LIBRARIES)
    add_subdirectory(${PBMM2_ThirdPartyDir}/minimap2 external/minimap2/build)
endif()
