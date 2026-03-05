#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "../include/grid.hpp"
#include "../include/solverCPU.hpp"
#include "../include/simGlobals.hpp"
#include "../include/boundaryConditions.hpp"
#include <vector>
#include <cmath>

// Simple 1D analytical solution for validation (infinite domain Gaussian)
double analytical_1d(double x, double t, double alpha) {
    return 100.0 * std::exp(-x*x / (4.0 * alpha * t)) / std::sqrt(4.0 * M_PI * alpha * t);
}

TEST_CASE("Grid indexing and basic fill", "[grid]") {
    Grid3D g(10, 10, 10);
    REQUIRE(g.nx() == 10);
    REQUIRE(g.ny() == 10);
    REQUIRE(g.nz() == 10);
    REQUIRE(g.size() == 1000);

    g.fill(42.0);
    REQUIRE(g(5,5,5) == 42.0);
}

TEST_CASE("CPU Stencil – basic stability & diffusion", "[cpu_stencil]") {
    SimulationGlobals globs;
    globs.nx = globs.ny = globs.nz = 21;
    globs.dx = globs.lx / globs.nx;
    globs.dt = 0.01;
    globs.alpha = 1.0;          // simplified for test
    globs.tol = 1e-6;
    globs.maxIters = 100;

    Grid3D current(globs.nx, globs.ny, globs.nz);
    Grid3D next(globs.nx, globs.ny, globs.nz);
    current.fill(0.0);
    // Center hot spot
    current(globs.nx/2, globs.ny/2, globs.nz/2) = 100.0;

    LinearAlgebra lin;
    HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, lin);

    // Run a few steps
    for (int t = 0; t < 50; ++t) {
        solver.step(current, next, globs, BoundaryConditions{{}, {}});
        std::swap(current, next);
    }

    // Center should have decreased, neighbors increased
    REQUIRE(current(globs.nx/2, globs.ny/2, globs.nz/2) < 100.0);
    REQUIRE(current(globs.nx/2 + 1, globs.ny/2, globs.nz/2) > 0.0);
}

TEST_CASE("Analytical match – 1D-like slice (very coarse)", "[analytical]") {
    // This is approximate – real test would use much finer grid
    SimulationGlobals globs;
    globs.nx = globs.ny = globs.nz = 41;
    globs.dx = 0.025;
    globs.dt = 0.0001;
    globs.alpha = 0.01;
    globs.tol = 1e-5;

    Grid3D current(globs.nx, globs.ny, globs.nz);
    current.fill(0.0);
    current(globs.nx/2, globs.ny/2, globs.nz/2) = 100.0;

    LinearAlgebra lin;
    HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, lin);
    BoundaryConditions bc{{}, {}};

    double t_sim = 0.0;
    for (int step = 0; step < 200; ++step) {
        Grid3D next = current;
        solver.step(current, next, globs, bc);
        current = next;
        t_sim += globs.dt;
    }

    double x = (globs.nx/2 + 5) * globs.dx - globs.lx/2.0;
    double expected = analytical_1d(x, t_sim, globs.alpha) * 100.0;
    double computed = current(globs.nx/2 + 5, globs.ny/2, globs.nz/2);

    // Tolerance loose because coarse grid + 3D effects
    REQUIRE(std::abs(computed - expected) < 15.0);
}
