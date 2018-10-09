#!/bin/bash
set -euo pipefail

# Main script
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash

module load gcc git ccache boost htslib ninja meson gtest zlib cram bedtools datamash samtools minimap2 gcovr pbbam pbcopper
meson -Db_coverage=true build .
ninja -C build test

cd build
find . -type f -iname '*.o' | xargs gcov -acbrfu {} \; >/dev/null && \
mkdir coverage && pushd coverage && mv ../*.gcov . && \
sed -i -e 's@Source:@Source:../@' *.gcov && \
sed -i -e 's@Graph:@Graph:../@' *.gcov && \
sed -i -e 's@Data:@Data:../@' *.gcov
