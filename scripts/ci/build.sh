#!/bin/bash
set -euo pipefail

# Main script
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash

module load gcc git ccache boost ninja meson seqan gtest pbcopper zlib htslib parasail cram bedtools datamash pbbam samtools
meson build-test .
ninja -C build-test

if [[ "${bamboo_planRepository_branchName}" == "develop" ]]; then
    echo "Copy binary to module"
    module purge
    module load gcc git ccache boost ninja meson seqan gtest htslib/static zlib/static
    LDFLAGS='-lrt' CXXFLAGS="-static-libgcc -static-libstdc++" meson --strip build-static . --prefix / -Dtests=false --default-library=static
    ninja -C build-static
    DESTDIR="/mnt/software/p/pbmm2/current" ninja -C build-static -v install
fi
