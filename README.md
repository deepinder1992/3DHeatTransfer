# 3DHeatTransfer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C++17-blue.svg)](https://isocpp.org/)
[![CUDA](https://img.shields.io/badge/CUDA-Enabled-76B900.svg)](https://developer.nvidia.com/cuda-toolkit)
[![CMake](https://img.shields.io/badge/CMake-3.20+-green.svg)](https://cmake.org/)

A high-performance **3D heat conduction solver** with four implicit solver backends (CPU and GPU). It supports **arbitrary complex geometries** using STL files and mixed Dirichlet/Neumann boundary conditions.

## Features

- Four implicit solver backends (all unconditionally stable):
  - **CPU Stencil** & **CUDA Stencil**: Stencil-based Jacobi iterations (lower memory usage)
  - **CPU Matrix** & **CUDA Matrix**: Explicit matrix + Jacobi solver (higher memory, more flexible)
- Arbitrary geometry support via STL files with automatic inlet/outlet/wall patch detection
- Mixed Dirichlet and Neumann boundary conditions
- VTK output for visualization in ParaView
- All command-line options have sensible defaults

## Quick Start

### 1. Build the project

```bash
git clone https://github.com/deepinder1992/3DHeatTransfer.git
cd 3DHeatTransfer
./build.sh

## 2. Run the solver
cd build

# Run with default settings
./heat3d

# Example commands
./heat3d --solver 1                    # CPU Stencil
./heat3d --solver 4 --nx 100           # CUDA Matrix with higher resolution
./heat3d --stlPath "../stlFiles/cube/cube.stl" # Stl file of choice

### 3. Command Line Options
	Option,Description,Default Value
	--solver,Solver type (1-4),3 (CUDA_STENCIL)
	--nx,Grid size (nx × nx × nx),50
	--steps,Maximum number of time steps,10000
	--dt,Time step size (seconds),1000
	--jacobiTol,Tolerance for Jacobi iterations,1e-6
	--globalTol,Global convergence tolerance,1e-8
	--verbosity,"Verbosity level (1=low, 2=med, 4=high)",1
	--writeInterval,Write VTK output every N steps,1000
	--blockDim,CUDA block size,216
	--bcTypeInlet,"Inlet BC (0=Dirichlet, 1=Neumann)",1
	--bcTypeOutlet,Outlet BC,1
	--bcTypeWall,Wall BC,0
	--bcValInlet,Inlet BC value,500.0
	--bcValOutlet,Outlet BC value,-500.0
	--bcValWall,Wall BC value,100.0
	--stlPath,Base STL geometry file path,../stlFiles/cylinder/cylinder.stl

### 4. Geometry (STL Files)
You must provide four STL files with the same base name:

name.stl
name_inlet.stl
name_outlet.stl
name_wall.stl

Example geometries are provided in the stlFiles/ folder.
All the files should be in one folder and you only need to provide the name of the "name.stl" file in the command 
rest will be picked automatically.

### 5. 
