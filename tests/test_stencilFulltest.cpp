#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/boundaryConditions.hpp"
#include "./include/tests_utils.hpp"
#include "../include/voxelReader.hpp"
#include<iostream>

Grid3D runCPUStencil() {

    Grid3D current(10,10,10,0.1);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CPU_STENCIL);

    Grid3D next = current;
    LinearAlgebra linAlgebra(50);
    BoundaryConditions bc(types, values);
    SimulationGlobals globs;

    HeatSolverCPUStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebra);

    for (int t = 0; t < 1000; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}

TEST(CPUStencilTest, FulltestStencilCPU ) {
    Grid3D current = runCPUStencil();

    EXPECT_NEAR(current(3,4,5), 100.587, 1e-2);
    EXPECT_NEAR(current(7,6,8), 109.268, 1e-2);
    EXPECT_NEAR(current(1,3,1), 91.795, 1e-2);

    EXPECT_NEAR(current(0,0,0), 98.100, 1e-2);
    EXPECT_NEAR(current(2,2,2), 94.300, 1e-2);
    EXPECT_NEAR(current(4,4,4), 97.789, 1e-2);
    EXPECT_NEAR(current(6,6,6), 103.023, 1e-2);
    EXPECT_NEAR(current(8,8,8), 104.624, 1e-2);

    EXPECT_NEAR(current(9,1,7), 100.0, 1e-2);
    EXPECT_NEAR(current(5,9,3), 100.0, 1e-2);
    EXPECT_NEAR(current(8,2,9), 114.193, 1e-2);

}

#ifdef ENABLE_CUDA
#include "../include/cudaHeaders/solverCUDA.cuh"
#include "../include/cudaHeaders/linearAlgebra.cuh"
#include "../include/cudaHeaders/boundaryConditions.cuh"

Grid3D runGPUStencil() {

    Grid3D current(10,10,10,0.1);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CUDA_STENCIL);

    Grid3D next = current;
    LinearAlgebraCUDA linAlgebraCUDA(50);
    BoundaryConditions bc(types, values);
    SimulationGlobals globs;

    HeatSolverCUDAStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebraCUDA);

    for (int t = 0; t < 1000; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}
TEST(CUDAStencilSolverTest, FulltestStencilGPU) {

    Grid3D current = runGPUStencil();

    EXPECT_NEAR(current(3,4,5), 100.587, 1e-2);
    EXPECT_NEAR(current(7,6,8), 109.268, 1e-2);
    EXPECT_NEAR(current(1,3,1), 91.795, 1e-2);

    EXPECT_NEAR(current(0,0,0), 98.100, 1e-2);
    EXPECT_NEAR(current(2,2,2), 94.300, 1e-2);
    EXPECT_NEAR(current(4,4,4), 97.789, 1e-2);
    EXPECT_NEAR(current(6,6,6), 103.023, 1e-2);
    EXPECT_NEAR(current(8,8,8), 104.624, 1e-2);

    EXPECT_NEAR(current(9,1,7), 100.0, 1e-2);
    EXPECT_NEAR(current(5,9,3), 100.0, 1e-2);
    EXPECT_NEAR(current(8,2,9), 114.193, 1e-2);

}


TEST(CompareCPUvsGPUStencilTest, CPUGPUEqualityStencil) {
     Grid3D gpuGridSten = runGPUStencil();
     Grid3D cpuGridSten = runCPUStencil();

    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < 10; ++k)
                EXPECT_NEAR(cpuGridSten(i, j, k), gpuGridSten(i, j, k), 1e-2);
}
#endif