#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/boundaryConditions.hpp"
#include "./include/tests_utils.hpp"
#include "../include/voxelReader.hpp"
#include<iostream>

Grid3D runCPUStencil() {

    Grid3D current(50,50,50);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CPU_STENCIL);

    Grid3D next = current;
    LinearAlgebra linAlgebra(50);
    
    SimulationGlobals globs;
    BoundaryConditions bc(globs.types, globs.values);

    HeatSolverCPUStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebra);

    for (int t = 0; t < 50; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}

TEST(CPUStencilTest, FulltestStencilCPU ) {
    Grid3D current = runCPUStencil();

    EXPECT_NEAR(current(45,45,45), 100.03799229143938, 1e-2);

    EXPECT_NEAR(current(10,5,0),   99.802401021107926, 1e-2);
    EXPECT_NEAR(current(20,10,5),  99.807292098812709, 1e-2);
    EXPECT_NEAR(current(30,15,10), 99.851957870753793, 1e-2);
    EXPECT_NEAR(current(40,20,15), 99.936802071288213, 1e-2);

    EXPECT_NEAR(current(45,25,20), 99.984087263715651, 1e-2);

    EXPECT_NEAR(current(0,25,10),  100.00000000000000, 1e-2);

    EXPECT_NEAR(current(5,30,20),  99.981302350706599, 1e-2);
    EXPECT_NEAR(current(10,35,25), 99.992330055712387, 1e-2);
    EXPECT_NEAR(current(15,40,30), 100.01235966601594, 1e-2);
    EXPECT_NEAR(current(20,45,35), 100.01959238164703, 1e-2);

}

#ifdef ENABLE_CUDA
#include "../include/cudaHeaders/solverCUDA.cuh"
#include "../include/cudaHeaders/linearAlgebra.cuh"
#include "../include/cudaHeaders/boundaryConditions.cuh"

Grid3D runGPUStencil() {

    Grid3D current(50,50,50);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CUDA_STENCIL);

    Grid3D next = current;
    LinearAlgebraCUDA linAlgebraCUDA(50);

    SimulationGlobals globs;
    BoundaryConditions bc(globs.types, globs.values);

    HeatSolverCUDAStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebraCUDA);

    for (int t = 0; t < 50; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}
TEST(CUDAStencilSolverTest, FulltestStencilGPU) {

    Grid3D current = runGPUStencil();

    EXPECT_NEAR(current(45,45,45), 100.03799229143938, 1e-2);

    EXPECT_NEAR(current(10,5,0),   99.802401021107926, 1e-2);
    EXPECT_NEAR(current(20,10,5),  99.807292098812709, 1e-2);
    EXPECT_NEAR(current(30,15,10), 99.851957870753793, 1e-2);
    EXPECT_NEAR(current(40,20,15), 99.936802071288213, 1e-2);

    EXPECT_NEAR(current(45,25,20), 99.984087263715651, 1e-2);

    EXPECT_NEAR(current(0,25,10),  100.00000000000000, 1e-2);

    EXPECT_NEAR(current(5,30,20),  99.981302350706599, 1e-2);
    EXPECT_NEAR(current(10,35,25), 99.992330055712387, 1e-2);
    EXPECT_NEAR(current(15,40,30), 100.01235966601594, 1e-2);
    EXPECT_NEAR(current(20,45,35), 100.01959238164703, 1e-2);
}


TEST(CompareCPUvsGPUStencilTest, CPUGPUEqualityStencil) {
     Grid3D gpuGridSten = runGPUStencil();
     Grid3D cpuGridSten = runCPUStencil();

    for (int i = 0; i < 50; ++i)
        for (int j = 0; j < 50; ++j)
            for (int k = 0; k < 50; ++k)
                EXPECT_NEAR(cpuGridSten(i, j, k), gpuGridSten(i, j, k), 1e-2);
}
#endif