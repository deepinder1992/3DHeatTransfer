# 3DHeatTransfer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C++17-blue.svg)](https://isocpp.org/)
[![CUDA](https://img.shields.io/badge/CUDA-Enabled-76B900.svg)](https://developer.nvidia.com/cuda-toolkit)
[![CMake](https://img.shields.io/badge/CMake-3.20+-green.svg)](https://cmake.org/)

## Summary

**3DHeatTransfer** is a high-performance solver for steady-state and transient heat conduction in three-dimensional domains. The software supports complex geometries defined via STL files and enables mixed Dirichlet and Neumann boundary conditions.

It provides both CPU and GPU implementations of implicit solvers, allowing users to balance flexibility, memory usage, and performance. The code is written in C++17 with optional CUDA acceleration and is designed for extensibility and reproducibility in computational heat transfer research.

---

## Statement of Need

Accurate simulation of heat conduction in complex geometries is essential in applications such as thermal management, energy systems, and materials processing. Existing general-purpose CFD tools are often heavyweight or not optimized for pure conduction problems.

**3DHeatTransfer** addresses this gap by providing:

* A lightweight, focused conduction solver
* Native support for STL-based geometries
* Multiple interchangeable solver backends (CPU/GPU, stencil/matrix)
* Straightforward integration into research workflows

---

## Features

* **Implicit solvers (unconditionally stable):**

  * CPU stencil-based Jacobi solver
  * CUDA stencil-based Jacobi solver
  * CPU matrix-based solver
  * CUDA matrix-based solver

* **Geometry handling:**

  * Arbitrary 3D geometries via STL files
  * Assign inlet, outlet, and wall patches

* **Boundary conditions:**

  * Mixed Dirichlet and Neumann conditions
  * Independent specification for inlet, outlet, and walls

* **Performance and usability:**

  * GPU acceleration using CUDA
  * Configurable solver parameters via command line
  * VTK output for visualization in ParaView

---

## Implementation

The solver uses a structured Cartesian grid and a finite-difference discretization of the heat equation.

Two approaches are implemented:

* **Stencil-based solvers:**
  Perform Jacobi iterations directly on the grid without assembling a global matrix, minimizing memory usage.

* **Matrix-based solvers:**
  Assemble the linear system explicitly and solve using Jacobi iterations, enabling greater flexibility at the cost of higher memory consumption.

CUDA implementations parallelize both approaches for execution on GPUs.

---

## Installation

### Requirements

* C++17 compatible compiler
* CMake ≥ 3.20
* CUDA toolkit (optional, for GPU solvers)
### Install Requirements
bash installDeps.sh   # Ubuntu
installDeps.bat       # Windows (MSYS2)
./build.sh

### Build

```bash
git clone https://github.com/deepinder1992/3DHeatTransfer.git
cd 3DHeatTransfer
./build.sh
```

---

## Usage

```bash
cd build
./heat3d
```

### Examples

```bash
# CPU stencil solver
./heat3d --solver 1

# CUDA matrix solver with higher resolution
./heat3d --solver 4 --nx 100

# Run with custom geometry
./heat3d --stlPath "../stlFiles/cube/cube.stl"
```

---

## Command Line Options

| Option            | Description                            | Default Value            |
| ----------------- | -------------------------------------- | ------------------------ |
| `--solver`        | Solver type (1–4)                      | 3 (CUDA_STENCIL)         |
| `--nx`            | Grid size (nx × nx × nx)               | 50                       |
| `--steps`         | Maximum number of time steps           | 10000                    |
| `--dt`            | Time step size (seconds)               | 1000                     |
| `--jacobiTol`     | Jacobi iteration tolerance             | 1e-6                     |
| `--globalTol`     | Global convergence tolerance           | 1e-8                     |
| `--verbosity`     | Verbosity level (1=low, 2=med, 4=high) | 1                        |
| `--writeInterval` | Write VTK output every N steps         | 1000                     |
| `--blockDim`      | CUDA block size                        | 216                      |
| `--bcTypeInlet`   | Inlet BC (0=Dirichlet, 1=Neumann)      | 1                        |
| `--bcTypeOutlet`  | Outlet BC                              | 1                        |
| `--bcTypeWall`    | Wall BC                                | 0                        |
| `--bcValInlet`    | Inlet BC value                         | 500.0                    |
| `--bcValOutlet`   | Outlet BC value                        | -500.0                   |
| `--bcValWall`     | Wall BC value                          | 100.0                    |
| `--stlPath`       | Base STL geometry file                 | ../stlFiles/cylinder/... |

---

## Geometry Input

The solver expects four STL files with a common base name:

```
name.stl
name_inlet.stl
name_outlet.stl
name_wall.stl
```

All files must be located in the same directory. Only the base file needs to be specified at runtime:

```bash
./heat3d --stlPath "../stlFiles/cylinder/cylinder.stl"
```

Associated boundary files are automatically detected.

---

## Output

Simulation results are written in VTK format and can be visualized using ParaView.

---
## Sample Output
Here the geometries simulated by the solver will be displayed

## Contributing

Contributions are welcome. Please open an issue or submit a pull request for improvements or bug fixes.

---

## License

This project is licensed under the MIT License.
