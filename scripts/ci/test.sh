#!/usr/bin/env bash
set -vex

########
# TEST #
########

ninja -C "${CURRENT_BUILD_DIR:-build}" -v test

############
# COVERAGE #
############

( cd "${CURRENT_BUILD_DIR:-build}" &&
  find . -type f -iname '*.o' | xargs gcov -acbrfu {} \; >/dev/null &&
    mkdir coverage && cd coverage && mv ../*.gcov . &&
    sed -i -e 's@Source:@Source:../@' *.gcov &&
    sed -i -e 's@Graph:@Graph:../@' *.gcov &&
    sed -i -e 's@Data:@Data:../@' *.gcov )
