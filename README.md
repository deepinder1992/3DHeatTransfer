# HeatTransfer3D

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C++17-blue.svg)](https://isocpp.org/)
[![CUDA](https://img.shields.io/badge/CUDA-Enabled-76B900.svg)](https://developer.nvidia.com/cuda-toolkit)
[![CMake](https://img.shields.io/badge/CMake-3.20+-green.svg)](https://cmake.org/)

---

## Summary

`HeatTransfer3D` is a lightweight, high-performance C++17 solver for steady-state and transient 3D heat conduction on structured Cartesian grids. 

The software offers a **streamlined and easy integration** from STL files to simulation by using a custom ray-tracing voxelization method that directly imports complex geometries and automatically classifies internal, external, and boundary regions. It supports mixed Dirichlet and Neumann boundary conditions on different patches.

A key feature is its **four interchangeable solver backends**: stencil-based and matrix-based solvers, each available on both CPU (with OpenMP) and CUDA (GPU). This multi-backend design allows users to easily trade off between memory usage, convergence speed, and hardware availability. Simulation results are exported in VTK format for easy visualization in ParaView.

HeatTransfer3D targets engineers and researchers in thermal management, electronics cooling, battery systems, additive manufacturing, and heat exchanger design who need fast, reproducible 3D heat conduction simulations with minimal setup overhead.


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
`HeatTransfer3D` uses a **factory pattern** to provide four solver backends, allowing users to choose the best balance between performance, convergence speed, and memory usage.

All solvers are based on the **Implicit Euler** time discretization scheme (unconditionally stable) and support both steady-state and transient simulations.

#### Solver Comparison

| Backend              | Time Scheme     | Linear Solver          | Solve Quality     | Memory Usage         | Best For                              |
|----------------------|-----------------|------------------------|-------------------|----------------------|---------------------------------------|
| **Stencil (Jacobi)** | Implicit Euler  | Jacobi iteration       | Approximate       | Very Low `O(1)`      | Large grids, GPU, memory-constrained runs |
| **Matrix (CG)**      | Implicit Euler  | Conjugate Gradient     | Near-exact        | Higher `O(N)`        | Faster convergence, smaller to medium problems |

#### Detailed Description

**Stencil-based Solvers (CPU + CUDA)**  
Apply a 7-point finite difference stencil directly on the temperature field and solve the implicit system using **Jacobi iteration**.  
- Extremely memory efficient  
- Highly optimized for GPU (coalesced memory access + shared memory)  
- Ideal for very large grids

**Matrix-based Solvers (CPU + CUDA)**  
Explicitly assemble a sparse coefficient matrix and solve the linear system using the **Conjugate Gradient** method.  
- Usually converges in fewer iterations than Jacobi  
- Higher memory usage (stores the matrix)  
- Better for cases where faster convergence is desired


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
## Project Structure

The main components are organized as follows:

- **`src/`**: Contains the core implementation
  - `main.cpp` — Command-line argument parsing and simulation workflow
  - `grid.cpp`, `voxelReader.cpp` — Grid management and STL voxelization
  - `boundaryConditions.cpp` — Boundary condition handling
  - `solverCPU.cpp`, `linearAlgebraCPU.cpp` — CPU solver implementations
  - **`cudaSrc/`** — All CUDA-specific code:
    - `solverCUDAStencil.cu`, `solverCUDAMatrix.cu`
    - `kernel.cu`, `linearAlgebraGPU.cu`, `boundaryConditions.cu`

- **`include/`**: Public headers and interfaces
  - Core abstractions: `solver.hpp`, `solverFactory.hpp`, `grid.hpp`, `voxelReader.hpp`
  - CPU-specific: `solverCPU.hpp`
  - CUDA-specific: `cudaHeaders/` directory containing `.cuh` files
  - Supporting classes: `sparseMatrix.hpp`, `heatMatrixBuilder.hpp`, `linearAlgebra.hpp`

- **`tests/`**: Comprehensive test suite that uses **GoogleTest (GTest)** and contains **30 tests** covering both CPU and GPU implementations. Key test files include:
  - `tests_analytical.cpp` — Analytical verification cases
  - `tests_grid.cpp`, `tests_VoxelReader.cpp` — Grid and voxelization tests
  - `tests_boundaryConditions.cpp` — Boundary condition validation
  - `tests_cpuLinearAlgebra.cpp`, `tests_cudaLinAlgebra.cu` — Linear algebra on both backends
  - `tests_sparseMatrix.cpp`
  - `test_stencilFulltest.cpp`, `test_matrixFulltest.cpp` — Full solver tests
  - `tests_vtkWriter.cpp` — Output writer validation
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

### Build

```bash
git clone https://github.com/deepinder1992/3DHeatTransfer.git
cd 3DHeatTransfer
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=OFF
cmake --build build --config Release -j4

```
To enable CUDA (if you have the CUDA Toolkit installed), use -DENABLE_CUDA=ON instead of -DENABLE_CUDA=OFF.

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
| `--nx`            | Grid size (nx × nx × nx)               | 100                       |
| `--steps`         | Maximum number of time steps           | 10000                    |
| `--dt`            | Time step size (seconds)               | 100                      |
| `--jacobiTol`     | Jacobi iteration tolerance             | 1e-6                     |
| `--globalTol`     | Global convergence tolerance           | 1e-8                     |
| `--verbosity`     | Verbosity level (1=low, 2=med, 4=high) | 1                        |
| `--writeInterval` | Write VTK output every N steps         | 1000                     |
| `--blockDim`      | CUDA block size                        | 512                      |
| `--conductivity`  | Conductivity of the material           | 10.0 W/m.K                |
| `--density`       | Density of the material                | 2200.0 kg/m3               |
| `--cp`            | Specific Heat of the material          | 800.0 J/kg.K               |
| `--bcTypeInlet`   | Inlet BC (0=Dirichlet, 1=Neumann)      | 0                        |
| `--bcTypeOutlet`  | Outlet BC                              | 0                        |
| `--bcTypeWall`    | Wall BC                                | 1                        |
| `--bcValInlet`    | Inlet BC value                         | 100.0                    |
| `--bcValOutlet`   | Outlet BC value                        | 100.0                   |
| `--bcValWall`     | Wall BC value                          | 100.0                    |
| `--stlPath`       | Base STL geometry file                 | ../stlFiles/cube/cube.stl |
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
![Cube](images/Cube.png)

**Figure 3**(Cube): Temperature distribution with Neumann boundary condition (heat flux) is applied on the four lateral faces, while Dirichlet conditions (fixed temperature) are imposed on the inlet and outlet faces located opposite to each other. (a) Cross-section parallel to the inlet–outlet direction. (b) Cross-section perpendicular to the inlet–outlet direction at the mid-plane.

![Cylinder](images/Cylinder.png)

**Figure 4**(Cylinder): Temperature field with Dirichlet conditions (100°C) at inlet and outlet, and Neumann heat flux on the curved wall. (a) Plane parallel to cylinder axis. (b) Plane perpendicular to the cylinder axis at the center height.

![SemiCylinder](images/SemiCylinder.png)

**Figure 6** (L-shaped channel): Temperature field with Dirichlet conditions imposed at the inlet on the upper horizontal face and at the outlet on the lower vertical face, while all remaining faces are treated as Neumann walls. (a) Plane through the L-shaped channel parallel to the primary flow path. (b) Plane perpendicular to the channel at center height.


## Contributing

Contributions are welcome. Please open an issue or submit a pull request for improvements or bug fixes.

---

## License

This project is licensed under the MIT License.
