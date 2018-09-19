#!/bin/bash
set -euo pipefail

# Main script
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash

module load gcc git ccache boost htslib ninja meson gtest zlib cram bedtools datamash samtools
meson build-test .
ninja -C build-test
ninja -C build-test test
