#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/boundaryConditions.hpp"
#include "./include/tests_utils.hpp"
#include "../include/voxelReader.hpp"
#include<iostream>


Grid3D runCPU_Matrix() {

    Grid3D current(50,50,50);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CPU_MATRIX);

    Grid3D next = current;

    LinearAlgebra linAlgebra(50);

    SimulationGlobals globs;
    BoundaryConditions bc(globs.types, globs.values);

    HeatSolverCPUMatrix solver(current, current.nx(), current.ny(), current.nz(),  globs.alpha, current.dx(),  globs.dt,  globs.k, bc, linAlgebra );

    for (int t = 0; t < 50; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}




TEST(CPUMatrixTest, FulltestMatrixCPU) {

    Grid3D current = runCPU_Matrix();

    EXPECT_NEAR(current(0,0,0),   99.989253181822718, 1e-2);
    EXPECT_NEAR(current(5,5,5),   99.939271760700052, 1e-2);
    EXPECT_NEAR(current(10,10,10),99.920903806242876, 1e-2);
    EXPECT_NEAR(current(20,20,20),99.960923805088171, 1e-2);

    EXPECT_NEAR(current(25,25,25),100.00453683973085, 1e-2);
    EXPECT_NEAR(current(30,30,30),100.04670867267922, 1e-2);
    EXPECT_NEAR(current(35,35,35),100.07408507124151, 1e-2);
    EXPECT_NEAR(current(40,40,40),100.07765875602330, 1e-2);
    EXPECT_NEAR(current(45,45,45),100.05372108215180, 1e-2);

    EXPECT_NEAR(current(10,5,0),   99.815877111838233, 1e-2);
    EXPECT_NEAR(current(20,10,5),  99.818735378592152, 1e-2);
    EXPECT_NEAR(current(30,15,10), 99.863842724720541, 1e-2);
    EXPECT_NEAR(current(40,20,15), 99.943372674442301, 1e-2);
    EXPECT_NEAR(current(45,25,20), 99.987207935252016, 1e-2);

    EXPECT_NEAR(current(0,25,10),  99.993601147456232, 1e-2);
    EXPECT_NEAR(current(5,30,20),  99.985260617931075, 1e-2);
    EXPECT_NEAR(current(10,35,25), 100.00241113424903, 1e-2);
    EXPECT_NEAR(current(15,40,30), 100.02666473797639, 1e-2);
    EXPECT_NEAR(current(20,45,35), 100.03331902866577, 1e-2);
}

#ifdef ENABLE_CUDA
#include "../include/cudaHeaders/solverCUDA.cuh"
#include "../include/cudaHeaders/linearAlgebra.cuh"
#include "../include/cudaHeaders/boundaryConditions.cuh"

Grid3D runGPU_Matrix() {

    Grid3D current(50,50,50);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CUDA_MATRIX);

    Grid3D next = current;

    LinearAlgebraCUDA linAlgebraCUDA(50);
    
    SimulationGlobals globs;
    BoundaryConditions bc(globs.types, globs.values);
    HeatSolverCUDAMatrix solver(current, current.nx(), current.ny(), current.nz(), globs.alpha, current.dx(), globs.dt, globs.k, bc, linAlgebraCUDA );

    for (int t = 0; t < 50; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}



TEST(CUDAMatrixTest, FulltestMatrixGPU) {

    Grid3D current = runGPU_Matrix();

    EXPECT_NEAR(current(0,0,0),   99.989253181822718, 1e-2);
    EXPECT_NEAR(current(5,5,5),   99.939271760700052, 1e-2);
    EXPECT_NEAR(current(10,10,10),99.920903806242876, 1e-2);
    EXPECT_NEAR(current(20,20,20),99.960923805088171, 1e-2);

    EXPECT_NEAR(current(25,25,25),100.00453683973085, 1e-2);
    EXPECT_NEAR(current(30,30,30),100.04670867267922, 1e-2);
    EXPECT_NEAR(current(35,35,35),100.07408507124151, 1e-2);
    EXPECT_NEAR(current(40,40,40),100.07765875602330, 1e-2);
    EXPECT_NEAR(current(45,45,45),100.05372108215180, 1e-2);

    EXPECT_NEAR(current(10,5,0),   99.815877111838233, 1e-2);
    EXPECT_NEAR(current(20,10,5),  99.818735378592152, 1e-2);
    EXPECT_NEAR(current(30,15,10), 99.863842724720541, 1e-2);
    EXPECT_NEAR(current(40,20,15), 99.943372674442301, 1e-2);
    EXPECT_NEAR(current(45,25,20), 99.987207935252016, 1e-2);

    EXPECT_NEAR(current(0,25,10),  99.993601147456232, 1e-2);
    EXPECT_NEAR(current(5,30,20),  99.985260617931075, 1e-2);
    EXPECT_NEAR(current(10,35,25), 100.00241113424903, 1e-2);
    EXPECT_NEAR(current(15,40,30), 100.02666473797639, 1e-2);
    EXPECT_NEAR(current(20,45,35), 100.03331902866577, 1e-2);
}

TEST(CompareCPUGPUMatrixTest, CPUGPUEqualityMatrix) {

    Grid3D cpu = runCPU_Matrix();
    Grid3D gpu = runGPU_Matrix();

    for (int i = 0; i < 50; ++i)
        for (int j = 0; j < 50; ++j)
            for (int k = 0; k < 50; ++k)
                EXPECT_NEAR(cpu(i,j,k), gpu(i,j,k), 1e-2);
}

#endif