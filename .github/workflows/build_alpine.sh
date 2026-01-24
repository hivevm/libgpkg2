#!/bin/sh

set -eu

cd "$WORKING_DIR"

apk add \
    g++ \
    make \
    cmake

mkdir build
cd build

cmake ..

make -j$(nproc)