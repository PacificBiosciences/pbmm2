#!/usr/bin/env bash
set -vex

###########
# INSTALL #
###########

if [[ ${PREFIX_ARG} ]]; then
  ## Cleaning out old installation from /mnt/software
  rm -rf "${PREFIX_ARG}"/*
fi

# Ensure code coverage and ASAN are disabled before installing
meson configure -Db_coverage=false -Db_sanitize=none "${CURRENT_BUILD_DIR:-build}"

DESTDIR="${DESTDIR:-/}" ninja -C "${CURRENT_BUILD_DIR:-build}" -v install
