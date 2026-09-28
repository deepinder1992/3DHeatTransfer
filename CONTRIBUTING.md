# Contributing to HeatTransfer3D

Thank you for your interest in contributing to HeatTransfer3D! Contributions of all kinds are welcome — bug reports, feature requests, documentation improvements, new test cases, and code contributions.

## How to Contribute

### Reporting Bugs

1. Check the [existing issues](https://github.com/deepinder1992/3DHeatTransfer/issues) to see if the bug has already been reported.
2. If not, open a new issue with:
   - A clear title and description
   - Steps to reproduce the problem
   - Expected vs. actual behaviour
   - Your platform (OS, compiler version, CUDA toolkit version if applicable)

### Suggesting Features

Open an issue describing:
- The problem the feature would solve
- A proposed approach (if you have one)
- Whether you are willing to implement it

### Submitting Code Changes

1. **Fork** the repository and create a feature branch from `main`:
   ```bash
   git checkout -b feature/my-improvement
   ```
2. **Build and test** your changes:
   ```bash
   cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug -DENABLE_CUDA=OFF
   cmake --build build
   cd build && ctest --output-on-failure
   ```
3. **Ensure all 30 existing tests pass** before submitting.
4. **Add tests** for any new functionality. Test files go in `tests/` and should use GoogleTest.
5. **Submit a pull request** against `main` with a clear description of what changed and why.

## Code Style

- C++17 standard
- Header files in `include/`, source files in `src/`, CUDA files in `src/cudaSrc/`
- Use descriptive variable and function names
- Keep functions focused — one responsibility per function
- Match the existing formatting style in the file you are editing

## Adding a New Solver Backend

1. Implement the `HeatSolver` interface defined in `include/solver.hpp` (provide a `step()` method).
2. Register the new backend in `include/solverFactory.hpp` with a unique integer flag.
3. Add corresponding tests in `tests/`.
4. Update the CLI help text in `src/main.cpp` and the README.

## Adding a New Boundary Condition Type

1. Extend the logic in `src/boundaryConditions.cpp` (for stencil solvers) and `include/heatMatrixBuilder.hpp` (for matrix solvers).
2. Add tests in `tests/tests_boundaryConditions.cpp`.

## Adding New Geometry

Simply place new STL files in `stlFiles/` following the naming convention:
```
name.stl          # Full geometry
name_inlet.stl    # Inlet patch
name_outlet.stl   # Outlet patch
name_wall.stl     # Wall patch
```

## Questions?

Open an issue on the [GitHub repository](https://github.com/deepinder1992/3DHeatTransfer/issues) or contact the maintainer directly.
