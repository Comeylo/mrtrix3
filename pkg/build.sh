#!/bin/bash#!/bin/bash

if [[ -d build ]]; then
    rm -rf build
fi
mkdir build
cd build

declare -a CMAKE_PLATFORM_FLAGS
if [[ ${HOST} =~ .*linux.* ]]; then
    CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/cross-linux.cmake")
fi

cmake .. \
      -DCMAKE_INSTALL_PREFIX="${PREFIX}" -DCMAKE_INSTALL_LIBDIR=lib  \
      ${CMAKE_PLATFORM_FLAGS[@]}

make -j${CPU_COUNT}
make install
