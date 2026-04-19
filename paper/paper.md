---
title: '3DHeatTransfer: A High-Performance Solver for 3D Heat Conduction with STL Geometry Support'
tags:
  - heat conduction
  - finite difference method
  - CUDA
  - GPU acceleration
  - STL geometry
  - Jacobi solver
authors:
  - name: Deepinder Singh
    orcid: **[EDIT: add your ORCID if you have one]**
    affiliation: 1
affiliations:
  - name: **[EDIT: Your institution / Independent / University of XYZ]**
    index: 1
date: **[EDIT: current date, e.g. 19 April 2026]**
bibliography: paper.bib
---

# Summary

**3DHeatTransfer** is a lightweight, high-performance C++ solver for steady-state and transient 3D heat conduction problems in complex geometries. It imports arbitrary 3D domains from STL files (with automatic patch detection for inlet/outlet/wall surfaces) and supports mixed Dirichlet and Neumann boundary conditions. 

Four interchangeable solver backends are provided: CPU and CUDA implementations of both stencil-based and matrix-based Jacobi iterative solvers. This design allows users to trade off memory usage and performance depending on available hardware. Results are exported in VTK format for easy visualization in ParaView.

The software targets researchers and engineers in thermal management, energy systems, materials processing, and related fields who need fast, focused conduction simulations without the overhead of general-purpose CFD packages.

# Statement of Need

Accurate modeling of heat conduction in complex three-dimensional geometries is critical in many scientific and engineering applications, including electronics cooling, additive manufacturing, battery thermal management, and heat exchanger design. While powerful general-purpose tools such as OpenFOAM or ANSYS exist, they can be heavyweight for pure conduction problems and often require extensive setup for simple thermal analyses.

**3DHeatTransfer** fills this gap by providing a specialized, easy-to-use, and performant open-source solver focused exclusively on the 3D heat equation. Its native STL geometry support and multiple solver backends (CPU/GPU, low-memory stencil vs. flexible matrix) make it particularly suitable for rapid parametric studies and integration into research workflows. The code emphasizes reproducibility, extensibility, and accessibility for users with varying hardware capabilities.

# State of the Field

Several open-source finite-difference and finite-volume heat transfer codes exist. However, many are either limited to simple Cartesian domains, lack GPU acceleration, or are embedded within much larger CFD frameworks. Tools with good STL support are often tied to complex meshing pipelines. **3DHeatTransfer** combines automatic STL patch detection, mixed boundary conditions, and four optimized solver variants in a single lightweight package.

# Software Design

The solver discretizes the heat equation on a uniform Cartesian grid using the finite-difference method. Two algorithmic approaches are implemented:

- **Stencil-based solvers**: Perform Jacobi iterations directly on the grid without assembling a global matrix (low memory footprint).
- **Matrix-based solvers**: Explicitly assemble and solve the sparse linear system (higher flexibility).

Both approaches have CPU and CUDA implementations, enabling seamless switching between backends via a command-line flag. Geometry is imported from four STL files (domain + separate inlet/outlet/wall patches). Boundary conditions and solver parameters (grid size, time step, tolerances, etc.) are fully configurable at runtime. VTK output ensures straightforward post-processing.

The code is written in modern C++17 with optional CUDA support and uses CMake for building. It has no heavy external dependencies beyond the CUDA toolkit when GPU support is enabled.

# Research Impact Statement

**3DHeatTransfer** enables faster iteration in thermal design studies by reducing the computational and setup overhead compared to general-purpose CFD software. Researchers can quickly evaluate different geometries, materials, and boundary conditions on both laptops (CPU) and high-performance GPU nodes. The open-source nature and MIT license facilitate community extensions, such as coupling with optimization routines or integration into larger multiphysics workflows.

# AI Usage Disclosure

# Acknowledgements

We thank the open-source community for the libraries and tools that made this project possible.

# References


