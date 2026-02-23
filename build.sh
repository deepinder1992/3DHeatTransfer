#!/bin/bash

BUILD_TYPE=${1:-Release}   # Default changed to Release

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

mkdir -p build
cd build || exit
cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE
cmake --build . -- -j$(nproc)
ls -l heat3d
