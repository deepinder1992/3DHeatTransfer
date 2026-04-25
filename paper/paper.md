---
title: 'HeatTransfer3D: A High-Performance Multi-Backend Solver for 3D Heat Conduction with STL Geometry Support'
tags:
  - heat conduction
  - finite difference method
  - CUDA
  - GPU acceleration
  - STL geometry
  - conjugate gradient
  - Implicit Jacobi
  - multi-backend
authors:
  - name: Deepinder Jot Singh Aulakh
    orcid: https://orcid.org/0000-0002-5140-2064
    affiliation: 1
affiliations:
  - name: Queen's University, Kingston, Canada
    index: 1
date: 20 April 2026
bibliography: paper.bib
---

# Summary

`HeatTransfer3D` is a lightweight, high-performance C++ solver for 3D heat conduction problems on structured Cartesian grids. It imports 3D geometries from STL files and incorporates them into the computational domain using voxelization. Geometry is mapped onto grid cells and used to define internal, external, and boundary regions. The solver supports mixed Dirichlet and Neumann boundary conditions on these regions.

The software provides **four interchangeable solver backends** — CPU and CUDA implementations of stencil-based and matrix-based Jacobi iterative solvers — allowing users to balance computational speed, memory usage, and hardware availability. Results are exported in VTK format for visualization in ParaView.

The software targets researchers and engineers working in thermal management, materials processing, electronics cooling, and related heat conduction applications.

# Statement of Need

Modeling three-dimensional heat conduction in complex geometries is a common requirement in applications such as electronics cooling, battery thermal management, additive manufacturing, and heat exchanger design. Despite its importance, existing tools often involve trade-offs between ease of use, computational efficiency, and geometric modeling complexity.

Commercial multiphysics platforms such as ANSYS and COMSOL Multiphysics provide comprehensive capabilities but are often inaccessible due to licensing costs and workflow complexity. Open-source frameworks such as OpenFOAM [@weller2007openfoam; @jasak2007openfoam] and FEniCS [@logg2012automated] offer high flexibility and general-purpose multiphysics functionality, but typically require substantial setup effort, including mesh generation and solver configuration, even for pure heat conduction problems.

In contrast, many lightweight finite-difference solvers are computationally efficient but are generally restricted to structured Cartesian grids and assume that the computational geometry is already represented on a grid. While GPU-accelerated implementations exist, they are typically focused on stencil optimization and performance improvement rather than providing an integrated workflow that connects complex geometries to structured-grid heat solvers.

`HeatTransfer3D` addresses this gap by providing a lightweight, focused, and high-performance open-source solver for steady-state and transient 3D heat conduction. The key feature of the software is the integration of geometry into a structured-grid framework using STL-based voxelization, where complex surfaces are mapped onto a Cartesian grid using ray tracing.

Key features include:
- Geometry handling via STL import with voxelization-based domain construction,
- Support for mixed Dirichlet and Neumann boundary conditions,
- Four interchangeable solver backends (CPU and CUDA implementations of stencil-based and matrix-based Jacobi solvers), enabling flexible trade-offs between performance, memory usage, and hardware availability,
- Simple command-line interface and VTK output for visualization in ParaView.

This makes the software suitable for rapid prototyping, parametric studies, and engineering applications where reproducibility, computational efficiency, and minimal setup overhead are critical.

# State of the Field
Existing tools for heat conduction simulation can be broadly categorized into three groups: general-purpose multiphysics frameworks, mesh-based finite element tools, and specialized research codes.

General-purpose frameworks such as OpenFOAM [@weller2007openfoam; @jasak2007openfoam] and FEniCS [@logg2012automated] provide extensive flexibility and support for coupled physics problems. However, they typically require substantial configuration effort, including mesh generation, solver selection, and case setup, which can be excessive for problems focused solely on heat conduction.

On the other end of the spectrum, lightweight finite-difference and finite-volume solvers are computationally efficient and widely used for heat conduction and diffusion-dominated problems due to their straightforward implementation on structured Cartesian grids [@leveque2007finite; @patankar1980numerical]. More advanced structured-grid frameworks, such as OpenSBLI [@howell2016opensbli], extend this class of methods to high-performance computing environments. However, such frameworks are designed as general-purpose PDE toolchains and typically require additional abstraction layers and configuration to define application-specific workflows.

Several GPU-accelerated finite-difference approaches have been developed to improve performance on structured grids. For example, phase-change heat conduction simulations on GPUs demonstrate substantial acceleration through optimized stencil execution, while remaining based on structured-grid discretizations that require complex geometries to be represented in grid-aligned form rather than directly incorporated from surface-based representations [@gpu_phasechange_heat]. Similarly, fast and interactive GPU-based heat conduction simulators have been proposed for two-dimensional problems, prioritizing real-time performance and interactivity over geometric modeling flexibility [@gpu_heat_2d_interactive]. These approaches are typically designed around specific research applications and optimized stencil implementations, and do not provide modular solver architectures for extensible simulation workflows.

In addition, simplified and educational solvers such as FDiff3 emphasize numerical understanding and clarity of implementation for heat conduction problems, primarily in pedagogical settings [@fdiff3]. Similarly, Python-based educational frameworks for two-dimensional heat transfer focus on accessibility and teaching purposes rather than extensible or modular solver design.

Meshless approaches such as radial basis function finite differences (RBF-FD) enable simulation on scattered nodes and can represent arbitrary three-dimensional geometries without structured grids [@fornberg2015solving; @miotti2021meshless]. However, these methods are typically CPU-based, involve increased formulation complexity, and require careful parameter selection.

HeatTransfer3D occupies a middle ground between these categories by combining:

- geometry-driven preprocessing that maps STL-based surfaces onto a structured Cartesian grid through voxelization,
- the simplicity and efficiency of finite-difference discretization on regular grids,
- a multi-backend design enabling both CPU and GPU execution,
- a lightweight implementation focused specifically on heat conduction problems.

This combination positions it as a focused engineering tool that integrates geometry-to-grid workflow handling with structured-grid solvers, rather than as a general-purpose multiphysics framework or an educational prototype.

# Software Design

`HeatTransfer3D` is implemented in C++17 with optional CUDA support for GPU acceleration. The codebase has a clean modular structure:

- `src/`: Main driver, geometry processing, boundary condition setup, solver kernels, and VTK output.
- `include/`: Core headers including `grid.hpp`, `boundaryConditions.hpp`, `sparseMatrix.hpp`, `heatMatrixBuilder.hpp`, `linearAlgebra.hpp`, `solver.hpp`, `solverCPU.hpp`, and `solverCUDA.hpp`.
- `tests/`: Unit and integration tests.
- `stlFiles/`: Example geometries with separate boundary patch files.

The steady 3D heat equation is discretized on a uniform Cartesian grid using the finite difference method. Geometry is incorporated from STL files via a ray-tracing voxelization routine. Each simulation requires four STL files: the main geometry plus three boundary patch files (`name_inlet.stl`, `name_outlet.stl`, `name_wall.stl`). Grid cells are classified as internal or assigned to one of the three boundary patches (inlet, outlet, wall).

**Boundary conditions** are applied as follows:
- Each patch (inlet, outlet, wall) can independently be set to Dirichlet (fixed temperature) or Neumann (fixed heat flux) using command-line flags `--bcTypeInlet`/`--bcTypeOutlet`/`--bcTypeWall` (0 = Dirichlet, 1 = Neumann) and the corresponding `--bcValXXX` values.
- These conditions are incorporated by modifying stencil coefficients (in stencil solvers) or the right-hand side and matrix entries (in matrix-based solvers) at boundary cells.

Four solver backends are available, selected at runtime with the `--solver` flag (1–4):

- **Stencil-based Jacobi solvers** (CPU with OpenMP, CUDA with coalesced memory access and shared memory): Perform Jacobi iterations directly on the temperature field using a 7-point stencil. These have a low memory footprint.
- **Matrix-based solvers** (CPU and CUDA): Explicitly assemble a sparse coefficient matrix using `sparseMatrix.hpp` and `heatMatrixBuilder.hpp`. The resulting linear system is solved using the **Conjugate Gradient** method implemented in `linearAlgebra.hpp` / `linearAlgebraCPU.cpp` (and the corresponding GPU version). This provides faster convergence compared to Jacobi iteration, especially for larger problems.

All solvers share a common interface defined in `solver.hpp`. Results are exported in legacy VTK format for visualization in ParaView.

The project builds with CMake. Helper scripts `build.sh` and `installDeps.sh` simplify compilation on Linux systems with optional CUDA support.

# Quality control

Correctness is ensured through unit tests, analytical verification cases, and convergence checks in the `tests/` directory. Tests validate STL voxelization accuracy, boundary patch labeling, correct application of mixed Dirichlet/Neumann conditions, and consistency of results across all four solver backends. Analytical solutions for 1D and 3D heat conduction problems are used to verify second-order spatial accuracy. The test suite can be executed with `ctest` after building.

# Performance

`HeatTransfer3D` offers four solver backends to balance speed and memory usage depending on hardware and problem size.

Performance was evaluated on a cube geometry with a residual tolerance of \(10^{-6}\).

![Execution time comparison](images/timing_bars.png)
**Figure 1:** Wall-clock time (s) for the four solvers on a \(100^3\) grid. CUDA backends show a clear advantage.

![Strong scaling](images/scaling_plot.png)

**Figure 2:** Strong scaling from \(50^3\) to \(150^3\) grid resolution.

Tests were performed on an Intel Core i5-11400H CPU and NVIDIA GeForce RTX 3050 GPU. The CUDA stencil backend delivers **8–20× speedup** over the CPU stencil version. Matrix-based Conjugate Gradient solvers provide faster convergence than Jacobi iterations for larger systems, while stencil solvers remain more memory-efficient.

# Simulations

The software has been tested on several sample geometries provided in the `stlFiles/` directory, including a cube, cylinder, L-shaped channel, and semi-cylinder. These cases demonstrate correct geometry import, boundary condition application, and solver convergence. Example temperature fields and convergence histories are shown in the repository documentation and can be reproduced using the supplied input files.

# Research Impact Statement

`HeatTransfer3D` fills an important gap by providing a lightweight, easy-to-use, yet high-performance tool for 3D steady-state heat conduction on complex STL geometries. Unlike full-featured CFD packages (e.g., OpenFOAM) that are heavy to install and learn, or simple 1D/2D educational codes that lack realistic geometry support, this software combines straightforward voxelization from separate inlet/outlet/wall STL files, flexible Dirichlet/Neumann boundary conditions per patch, and multiple high-performance solvers (stencil Jacobi and matrix-based Conjugate Gradient on both CPU and GPU) in a single package.

The ability to switch between low-memory stencil solvers and faster-converging Conjugate Gradient matrix solvers enables users to balance speed and memory usage depending on hardware and problem size. This makes the software particularly valuable for:

- Parametric thermal studies in electronics cooling and heat sink design
- Battery thermal management simulations
- Materials processing and casting applications
- Rapid prototyping of thermal designs before moving to more expensive full-physics solvers

This solver lowers the barrier for researchers and engineers to perform reproducible 3D heat transfer simulations. The modular design (especially the shared sparse matrix infrastructure) also makes it straightforward to extend with additional physics modules in the future. By releasing the software openly, we hope to encourage community contributions and adoption in both academic research and industrial workflows where fast turnaround on moderately complex geometries is needed.

# AI Usage Disclosure

AI-based tools were used to assist with debugging, syntax support, and language refinement during the development and writing of this work. All scientific decisions, implementation design, and results validation were performed and verified by the authors.

The following tools were used:
- ChatGPT (GPT-4o)
- Grok
  
# Availability

The source code is openly available under the MIT license at  
https://github.com/deepinder1992/3DHeatTransfer.  

Installation is performed via CMake on Linux systems supporting C++17 and optional CUDA. Detailed instructions, including dependency installation and build commands, are provided in the repository README.

# Acknowledgements

The author thanks the open-source community for the tools that made this project possible.

# References
