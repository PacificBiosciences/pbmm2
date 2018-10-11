#!/bin/bash
set -euo pipefail

# Main script
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash

module load gcc git ccache boost htslib ninja meson gtest zlib cram bedtools datamash samtools minimap2
meson build-test .
ninja -C build-test test
