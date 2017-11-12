##########################################
# CMake build script for the PBMM2 aligner
##########################################

cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0048 NEW)
project(UNANIMITY VERSION 0.0.1 LANGUAGES CXX C)

set(ROOT_PROJECT_NAME ${PROJECT_NAME} CACHE STRING "root project name")

# Build type
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile RelWithDebInfo ReleaseWithAssert" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Main project paths
set(PBMM2_RootDir       ${UNANIMITY_SOURCE_DIR})
set(PBMM2_ThirdPartyDir ${PBMM2_RootDir}/third-party)

# Project configuration
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

# Fixed order, do not sort or shuffle
include(pbmm2-releasewithassert)
include(pbmm2-dependencies)
include(pbmm2-compilerflags)

include(${PBMM2_RootDir}/pbmm2.cmake)
