#!/usr/bin/env bash

set -e

echo "====================================="
echo "Installing dependencies for Heat3D"
echo "Ubuntu setup"
echo "====================================="

sudo apt-get update

sudo apt-get install -y \
    build-essential \
    cmake \
    git \
    libomp-dev \
    libcli11-dev

echo ""
echo "NOTE:"
echo "- CUDA must be installed separately from NVIDIA toolkit"
echo "- Ensure nvcc is in PATH before building GPU targets"
echo ""

echo "Done."
echo "Next: ./build.sh"
