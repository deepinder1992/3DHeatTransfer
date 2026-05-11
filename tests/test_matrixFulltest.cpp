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

    EXPECT_NEAR(current(0,0,0),     99.430742224719765, 1e-2);

    EXPECT_NEAR(current(5,5,5),     83.331495453146303, 1e-2);
    EXPECT_NEAR(current(10,10,10),  69.659287900184253, 1e-2);
    EXPECT_NEAR(current(20,20,20),  55.175265483002114, 1e-2);

    EXPECT_NEAR(current(25,25,25),  53.820026070704820, 1e-2);
    EXPECT_NEAR(current(30,30,30),  55.869669403381678, 1e-2);
    EXPECT_NEAR(current(35,35,35),  61.756596929444378, 1e-2);
    EXPECT_NEAR(current(40,40,40),  72.085561261875583, 1e-2);
    EXPECT_NEAR(current(45,45,45),  86.438196025453664, 1e-2);

    EXPECT_NEAR(current(10,5,0),    98.408678165825179, 1e-2);
    EXPECT_NEAR(current(20,10,5),   82.416451471664729, 1e-2);
    EXPECT_NEAR(current(30,15,10),  69.090936738369109, 1e-2);
    EXPECT_NEAR(current(40,20,15),  60.468221796308896, 1e-2);
    EXPECT_NEAR(current(45,25,20),  56.224065882411807, 1e-2);

    EXPECT_NEAR(current(0,25,10),   70.565215606578846, 1e-2);

    EXPECT_NEAR(current(5,30,20),   56.133293526213478, 1e-2);
    EXPECT_NEAR(current(10,35,25),  54.482955977004679, 1e-2);
    EXPECT_NEAR(current(15,40,30),  56.464368402598048, 1e-2);
    EXPECT_NEAR(current(20,45,35),  62.468401137723113, 1e-2);
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

    EXPECT_NEAR(current(0,0,0),     99.430742224719765, 1e-2);

    EXPECT_NEAR(current(5,5,5),     83.331495453146303, 1e-2);
    EXPECT_NEAR(current(10,10,10),  69.659287900184253, 1e-2);
    EXPECT_NEAR(current(20,20,20),  55.175265483002114, 1e-2);

    EXPECT_NEAR(current(25,25,25),  53.820026070704820, 1e-2);
    EXPECT_NEAR(current(30,30,30),  55.869669403381678, 1e-2);
    EXPECT_NEAR(current(35,35,35),  61.756596929444378, 1e-2);
    EXPECT_NEAR(current(40,40,40),  72.085561261875583, 1e-2);
    EXPECT_NEAR(current(45,45,45),  86.438196025453664, 1e-2);

    EXPECT_NEAR(current(10,5,0),    98.408678165825179, 1e-2);
    EXPECT_NEAR(current(20,10,5),   82.416451471664729, 1e-2);
    EXPECT_NEAR(current(30,15,10),  69.090936738369109, 1e-2);
    EXPECT_NEAR(current(40,20,15),  60.468221796308896, 1e-2);
    EXPECT_NEAR(current(45,25,20),  56.224065882411807, 1e-2);

    EXPECT_NEAR(current(0,25,10),   70.565215606578846, 1e-2);

    EXPECT_NEAR(current(5,30,20),   56.133293526213478, 1e-2);
    EXPECT_NEAR(current(10,35,25),  54.482955977004679, 1e-2);
    EXPECT_NEAR(current(15,40,30),  56.464368402598048, 1e-2);
    EXPECT_NEAR(current(20,45,35),  62.468401137723113, 1e-2);
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