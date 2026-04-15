# 3DHeatTransfer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C++17-blue.svg)](https://isocpp.org/)
[![CUDA](https://img.shields.io/badge/CUDA-Enabled-76B900.svg)](https://developer.nvidia.com/cuda-toolkit)
[![CMake](https://img.shields.io/badge/CMake-3.20+-green.svg)](https://cmake.org/)

A high-performance **3D heat conduction solver** supporting CPU and GPU (CUDA) with both explicit (stencil) and implicit (matrix + Jacobi) methods. It supports **arbitrary complex geometries** using STL files and mixed Dirichlet/Neumann boundary conditions.

## Features

- Four solver backends:
  - 1 → CPU Stencil (explicit)
  - 2 → CPU Matrix (implicit)
  - 3 → CUDA Stencil (explicit)
  - 4 → CUDA Matrix (implicit)
- Arbitrary geometry support via STL files with automatic inlet/outlet/wall patch detection
- Mixed boundary conditions (Dirichlet and Neumann)
- VTK output for visualization in ParaView
- All options have sensible defaults

## Quick Start

### Build the project

```bash
git clone https://github.com/deepinder1992/3DHeatTransfer.git
cd 3DHeatTransfer
./build.sh
