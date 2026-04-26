#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include "tests_utils.hpp"
#include "../include/grid.hpp"
#include "../include/simGlobals.hpp"
#include "../include/boundaryConditions.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"

// Helper: 1D Gaussian analytical solution (infinite domain)
double analytical_1d(double x, double t, double alpha) {
    if (t <= 0.0) return (std::abs(x) < 1e-10) ? 100.0 : 0.0;
    return 100.0 * std::exp(-x*x / (4.0 * alpha * t)) / std::sqrt(4.0 * M_PI * alpha * t);
}

bool approx_equal(double a, double b, double tol = 1e-4) {
    return std::abs(a - b) < tol;
}

#ifdef ENABLE_CUDA
int main_cuda_tests();
#endif

int main() {
    std::cout << "=== Running HeatTransfer3D Tests ===\n\n";

    // Test 1: Grid basics
    {
        std::cout << "Test 1: Grid indexing & fill... ";
        Grid3D g(10, 12, 8, 0.1);
        if (g.nx() != 10 || g.ny() != 12 || g.nz() != 8 || g.size() != 960) {
            std::cout << "FAIL\n";
            return 1;
        }
        g.fill(42.5);
        if (!approx_equal(g(5,6,4), 42.5)) {
            std::cout << "FAIL\n";
            return 1;
        }
        std::cout << "PASS\n";
    }

    // Test 2: CPU stencil diffusion (qualitative)
    {
        std::cout << "Test 2: Stencil diffusion (qualitative)... ";
        SimulationGlobals globs;
        globs.nx = globs.ny = globs.nz = 21;
        globs.dx = 1.0 / globs.nx;
        globs.dt = 0.005;
        globs.alpha = 0.1;
        globs.tol = 1e-5;
        globs.maxIters = 200;

        Grid3D current(globs.nx, globs.ny, globs.nz, globs.dx);
        current.fill(0.0);
        current(globs.nx/2, globs.ny/2, globs.nz/2) = 100.0;

        LinearAlgebra lin(globs.maxIters);
        HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, lin);
        BoundaryConditions bc{{}, {}};

        for (int step = 0; step < 50; ++step) {
            Grid3D next = current;
            solver.step(current, next, globs, bc);
            current = next;
        }

        double center = current(globs.nx/2, globs.ny/2, globs.nz/2);
        double neighbor = current(globs.nx/2 + 1, globs.ny/2, globs.nz/2);

        if (center >= 100.0 || neighbor <= 0.0) {
            std::cout << "FAIL (center=" << center << ", neighbor=" << neighbor << ")\n";
            return 1;
        }
        std::cout << "PASS (center cooled to " << center << ")\n";
    }

    // Test 3: Rough analytical match
    {
        std::cout << "Test 3: Analytical approximation... ";
        SimulationGlobals globs;
        globs.nx = globs.ny = globs.nz = 41;
        globs.dx = 0.025;
        globs.dt = 0.00005;
        globs.alpha = 0.01;
        globs.tol = 1e-4;
        globs.maxIters = 300;

        Grid3D current(globs.nx, globs.ny, globs.nz, globs.dx);
        current.fill(0.0);
        current(globs.nx/2, globs.ny/2, globs.nz/2) = 100.0;

        LinearAlgebra lin(globs.maxIters);
        HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, lin);
        BoundaryConditions bc{{}, {}};

        double t_sim = 0.0;
        for (int step = 0; step < 400; ++step) {
            Grid3D next = current;
            solver.step(current, next, globs, bc);
            current = next;
            t_sim += globs.dt;
        }

        double x = 5.0 * globs.dx;
        double expected = analytical_1d(x, t_sim, globs.alpha);
        double computed = current(globs.nx/2 + 5, globs.ny/2, globs.nz/2);

        if (std::abs(computed - expected) > 20.0) {
            std::cout << "FAIL (computed=" << computed << ", expected≈" << expected << ")\n";
            return 1;
        }
        std::cout << "PASS (difference within loose tolerance)\n";
    }

    std::cout << "All CPU tests passed!\n";

    #ifdef ENABLE_CUDA
        main_cuda_tests();
    #endif
    return 0;
}
