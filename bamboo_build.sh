#!/usr/bin/env bash
set -e

################
# DEPENDENCIES #
################

## Load modules
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

module load gcc
module load ccache
module load git

module load meson
module load ninja

module load boost
module load cram
module load gtest
module load minimap2
module load bedtools
module load samtools
module load datamash
module load gcovr

case "${bamboo_planRepository_branchName}" in
  master)
    module load pbbam/master
    module load pbcopper/master
    ;;
  *)
    module load pbbam/develop
    module load pbcopper/develop
    ;;
esac
set -vx

BOOST_ROOT="${BOOST_ROOT%/include}"
# unset these variables to have meson discover all
# boost-dependent variables from BOOST_ROOT alone
unset BOOST_INCLUDEDIR
unset BOOST_LIBRARYDIR

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ $USER == bamboo ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi

case "${bamboo_planRepository_branchName}" in
  develop|master)
    export PREFIX_ARG="/mnt/software/p/pbmm2/${bamboo_planRepository_branchName}"
    export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
    ;;
  *)
    export BUILD_NUMBER="0"
    ;;
esac

# call the main build+test scripts
export ENABLED_TESTS="true"
export LDFLAGS="-static-libstdc++ -static-libgcc"

source scripts/ci/build.sh
source scripts/ci/test.sh

if [[ -z ${PREFIX_ARG+x} ]]; then
  echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
  exit 0
fi

source scripts/ci/install.sh
