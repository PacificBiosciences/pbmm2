#!/usr/bin/env bash
set -vex

export BUILD_NUMBER="0"
export ENABLED_TESTS="true"

case "${GCC_VERSION}" in
  4.8)
    module load gtest/gcc48

    # have to use wrap fallback due to ABI changes
    module unload pbbam pbcopper
    export WRAP_MODE="forcefallback"

    module load zlib
    module load htslib
    ;;

  next)
    module load gcc/8.1.0
    module load gtest
    ;;

  *)
    case "${bamboo_planRepository_branchName}-${bamboo_shortJobKey}" in
      develop-GM|master-GM)
        export PREFIX_ARG="/mnt/software/p/pbmm2/${bamboo_planRepository_branchName}"
        export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
        ;;
    esac

    module load gcc
    module load gtest
    ;;
esac

module load ccache

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ $USER == bamboo ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.${bamboo_shortJobKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
