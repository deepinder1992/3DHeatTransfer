# HeatTransfer3D

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C++17-blue.svg)](https://isocpp.org/)
[![CUDA](https://img.shields.io/badge/CUDA-Enabled-76B900.svg)](https://developer.nvidia.com/cuda-toolkit)
[![CMake](https://img.shields.io/badge/CMake-3.20+-green.svg)](https://cmake.org/)

---

## Summary

`HeatTransfer3D` is a high-performance solver for steady-state and transient heat conduction in three-dimensional domains on structured Cartesian grids. The software imports geometries defined via STL files and maps them onto the computational grid using a voxelization procedure, where grid cells are classified into internal, external, and boundary regions.

The solver supports mixed Dirichlet and Neumann boundary conditions and provides both CPU and GPU implementations of iterative solvers. It is written in C++17 with optional CUDA acceleration and is designed for reproducibility and extensibility in heat conduction simulations.

---

## Statement of Need

Accurate simulation of heat conduction in complex geometries is essential in applications such as thermal management, energy systems, and materials processing. Existing general-purpose CFD tools are often complex and not optimized for standalone heat conduction problems.

**HeatTransfer3D** addresses this by providing:
- A lightweight structured-grid heat conduction solver
- STL-based geometry handling via voxelization onto Cartesian grids
- Mixed boundary condition support (Dirichlet and Neumann)
- Multiple interchangeable solver backends (CPU/GPU, stencil-based and matrix-based Jacobi iterations)

This design enables efficient simulation workflows while maintaining a simple and reproducible setup for research and engineering applications.

---

## Features

### Solvers
- CPU stencil-based Jacobi iteration
- CUDA stencil-based Jacobi iteration
- CPU matrix-based Jacobi solver
- CUDA matrix-based Jacobi solver

### Geometry Handling
- STL-based geometry input
- Geometry mapped onto structured grids using voxelization
- Boundary regions defined as inlet, outlet, and wall through grid classification

### Boundary Conditions
- Mixed Dirichlet and Neumann conditions
- Independent specification for inlet, outlet, and wall regions

### Performance and Usability
- GPU acceleration via CUDA
- Configurable solver parameters via command line
- VTK output for visualization in ParaView

---

## Implementation

The solver is based on a structured Cartesian grid discretization of the heat equation.

Two iterative approaches are implemented:

- **Stencil-based solvers**  
  Jacobi iterations are performed directly on the grid without explicit matrix assembly, minimizing memory overhead.

- **Matrix-based solvers**  
  The discretized system is explicitly assembled and solved using Jacobi iterations, enabling more explicit control at the cost of higher memory usage.

CUDA implementations accelerate both approaches on GPUs.

---

## Installation

### Requirements

* C++17 compatible compiler
* CMake ≥ 3.20
* CUDA toolkit (optional, for GPU solvers)
### Install Requirements
Run the following commands to install the requirements for the code to build and run
  1. **Ubuntu**:    bash installDeps.sh   
  2. **Windows (MSYS2)**: installDeps.bat       

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
|  `-h,--help`      | Print this summary                     | NA                       |

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
## Visualizing Results in ParaView

Simulation results are exported as `.vti` (VTK Image Data) files on a full Cartesian grid.

To view the actual geometry in ParaView:

- Open the `.vti` file  
- Apply a **Threshold** filter  
- Select the appropriate mask field (e.g., domain/solid indicator) to extract the region of interest  

After filtering, you can visualize temperature or other fields as usual.

---
## Sample Output
Here the geometries simulated by the solver will be displayed

## Contributing

Contributions are welcome. Please open an issue or submit a pull request for improvements or bug fixes.

---

## License

This project is licensed under the MIT License.
