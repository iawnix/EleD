#!/usr/bin/env bash
set -euo pipefail

# Build xTB as a shared library for the optional ctypes backend.
# Usage: build_xtb_shared.sh /path/to/xtb-source /path/to/prefix

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 XTB_SOURCE_DIR [INSTALL_PREFIX]" >&2
  exit 2
fi

source_dir=$(realpath "$1")
prefix=$(realpath -m "${2:-$PWD/.xtb}")
build_dir="$prefix/build"

if [[ ! -f "$source_dir/CMakeLists.txt" ]]; then
  echo "xTB source directory does not contain CMakeLists.txt: $source_dir" >&2
  exit 2
fi

cmake -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build "$build_dir" --parallel
cmake --install "$build_dir"

mapfile -t shared_libraries < <(find "$prefix" -type f -name 'libxtb.so*' -print)
if [[ ${#shared_libraries[@]} -eq 0 ]]; then
  echo "xTB build completed without installing libxtb.so" >&2
  echo "Inspect the xTB CMake options and enable its shared-library target." >&2
  exit 1
fi

echo "Installed xTB shared libraries:"
printf '  %s\n' "${shared_libraries[@]}"
echo "Use XTB_LIBRARY=${shared_libraries[0]} with eled's native backend."
