---
title: '3DHeatTransfer: A High-Performance Multi-Backend Solver for 3D Heat Conduction with STL Geometry Support'
tags:
  - heat conduction
  - finite difference method
  - CUDA
  - GPU acceleration
  - STL geometry
  - Jacobi solver
  - multi-backend
authors:
  - name: Deepinder Singh
    orcid: **[EDIT: add your ORCID if available]**
    affiliation: 1
affiliations:
  - name: **[EDIT: Your affiliation / institution]**
    index: 1
date: 19 April 2026
bibliography: paper.bib
---

# Summary

**3DHeatTransfer** is a lightweight, high-performance C++ solver for steady-state and transient 3D heat conduction problems in complex geometries. It imports arbitrary 3D domains from STL files (with patch assignment: inlet, outlet, and  wall) and supports mixed Dirichlet and Neumann boundary conditions. 

The software provides **four interchangeable solver backends** — CPU and CUDA implementations of both stencil-based and matrix-based Jacobi iterative solvers — allowing users to select the best combination of speed, memory usage, and hardware availability. Results are exported in VTK format for easy visualization in ParaView.

The software targets researchers and engineers working in thermal management, energy systems, materials processing, electronics cooling, and related fields.

# Statement of Need

Accurate modeling of heat conduction in complex three-dimensional geometries is critical in many scientific and engineering applications, including electronics cooling, additive manufacturing, battery thermal management, and heat exchanger design. While powerful general-purpose tools such as OpenFOAM or ANSYS exist, they can be heavyweight for pure conduction problems.

**3DHeatTransfer** fills this gap by providing a specialized, easy-to-use, and highly performant open-source solver focused exclusively on the 3D heat equation, with flexible multi-backend solver selection and native STL geometry support.

# Software Design

The codebase follows a clean, modular structure written in modern C++17 with optional CUDA support. The main directories are organized as follows:

- **`src/`**: Core implementation files, including the main driver, geometry processing, boundary condition handling, solver kernels, and VTK output routines.
- **`include/`**: Header files defining classes and functions for the grid, solvers, STL importer, and utilities.
- **`tests/`**: Unit and integration tests.
- **`stlFiles/`**: Sample geometries (cube, cylinder, L_Channel, semiCylinder) with associated boundary patch files.

The solver uses a uniform Cartesian grid with finite-difference discretization of the heat equation. Two fundamental algorithmic approaches are implemented:

1. **Stencil-based solvers** — Direct Jacobi iterations on the temperature field (low memory footprint, high performance).
2. **Matrix-based solvers** — Explicit assembly of a sparse coefficient matrix followed by Jacobi iterations (more flexible for extensions).

Each approach has both a CPU and a highly optimized CUDA (GPU) implementation. Users can select any of the four backends at runtime using a simple command-line flag (`--solver 1..4`). CUDA kernels are optimized for coalesced memory access and make use of shared memory where beneficial. The geometry pipeline automatically detects and labels boundary patches from separate STL files, simplifying the application of different boundary conditions.

CMake is used for building, with convenient `build.sh` and `installDeps` scripts provided.

# Performance

One of the key strengths of **3DHeatTransfer** is its **multi-backend design**, which lets users dynamically choose the optimal solver for their hardware and problem size.

Performance comparisons on a representative test case (cube geometry, steady-state convergence) are shown below:

![CPU vs GPU Speedup](performance/cpu_vs_gpu_speedup.png)  
**Figure 1:** Execution time comparison of the four solvers for a 100×100×100 grid (lower is better). CUDA backends demonstrate substantial speedup over CPU versions.

![Strong Scaling](performance/strong_scaling.png)  
**Figure 2:** Strong scaling with increasing grid resolution (50³ to 200³). The GPU stencil backend maintains excellent performance thanks to its low memory overhead.

On typical hardware, the CUDA stencil solver achieves **8–25× speedup** compared to the CPU stencil solver for grids larger than 80³, while still offering good performance on laptops via the CPU backends. This flexibility makes the software suitable for both rapid prototyping and large-scale simulations.

# State of the Field

Several open-source heat transfer codes exist, but most are either limited to simple domains, lack GPU acceleration, or are part of much larger CFD frameworks. **3DHeatTransfer** stands out by combining automatic STL support, mixed boundary conditions, and four optimized solver backends in a single lightweight and easy-to-use package.

# Research Impact Statement

By offering high performance together with flexible solver selection, **3DHeatTransfer** enables faster design iterations and parametric studies in thermal engineering research. The open-source MIT license and modular structure also facilitate community contributions and integration into larger multiphysics workflows.

# AI Usage Disclosure

No AI tools were used in the development or writing of this software/paper beyond standard grammar and spell-checking assistance.

# Acknowledgements

We thank the open-source community for the tools that made this project possible.

# References
