#!/usr/bin/env bash
set -e

################
# DEPENDENCIES #
################

## Load modules
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

module load git

module load meson
module load ninja

module load boost
module load cram

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

# call the main build+test scripts
export LDFLAGS="-static-libstdc++ -static-libgcc"

source scripts/ci/setup.sh
source scripts/ci/build.sh
source scripts/ci/test.sh

if [[ ${BUILD_NUMBER} == 0 ]]; then
  echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
  exit 0
fi

source scripts/ci/install.sh
