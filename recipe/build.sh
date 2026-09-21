#!/usr/bin/env bash
set -euxo pipefail

cmake -S "$SRC_DIR" -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"

cmake --build build --parallel "${CPU_COUNT:-4}"
cmake --install build
