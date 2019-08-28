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

if [[ -z ${bamboo_planRepository_branchName+x} ]]; then
  : #pass
elif [[ ! -d /pbi/flash/bamboo/ccachedir ]]; then
  echo "[WARNING] /pbi/flash/bamboo/ccachedir is missing"
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.develop
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $bamboo_planRepository_branchName == master ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.master
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $USER == bamboo ]]; then
  _shortPlanKey=$(echo ${bamboo_shortPlanKey}|sed -e 's/[0-9]*$//')
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}
  if [[ -d /pbi/flash/bamboo/ccachedir/${_shortPlanKey}.${bamboo_shortJobKey}.develop ]]; then
    cp -a /pbi/flash/bamboo/ccachedir/${_shortPlanKey}.${bamboo_shortJobKey}.develop $CCACHE_DIR
  fi
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
