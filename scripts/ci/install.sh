#!/usr/bin/env bash
set -vex

###########
# INSTALL #
###########

if [[ ${PREFIX_ARG} ]]; then
  ## Cleaning out old installation from /mnt/software
  rm -rf "${PREFIX_ARG}"/*
fi

# Ensure code coverage is disabled before installing
meson configure -Db_coverage=false "${CURRENT_BUILD_DIR:-build}"

DESTDIR="${DESTDIR:-/}" ninja -C "${CURRENT_BUILD_DIR:-build}" -v install
