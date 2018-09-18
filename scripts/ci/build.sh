#!/bin/bash
set -euo pipefail

# Main script
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash

module load gcc git ccache boost htslib/1.8 ninja meson gtest zlib cram bedtools datamash
meson build-test .
ninja -C build-test

if [[ "${bamboo_planRepository_branchName}" == "develop" ]]; then
    echo "Copy binary to module"
    module purge
    module load htslib/1.8-static zlib/static gcc git ccache boost ninja meson gtest
    # LDFLAGS='-lrt' CXXFLAGS="-static-libgcc -static-libstdc++" meson --strip build-static . --prefix / -Dtests=false --default-library=static
    meson --strip build-static . --prefix / -Dtests=false
    ninja -C build-static
    DESTDIR="/mnt/software/p/pbmm2/current" ninja -C build-static -v install
fi
