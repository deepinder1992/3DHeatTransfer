#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/boundaryConditions.hpp"
#include "./include/tests_utils.hpp"
#include "../include/voxelReader.hpp"
#include<iostream>


Grid3D runCPU_Matrix() {

    Grid3D current(10,10,10,0.1);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CPU_MATRIX);

    Grid3D next = current;

    LinearAlgebra linAlgebra(50);
    BoundaryConditions bc(types, values);
    SimulationGlobals globs;

    HeatSolverCPUMatrix solver(
        current,
        current.nx(),
        current.ny(),
        current.nz(),
        1e-4,
        current.dx(),
        globs.dt,
        globs.k,
        bc,
        linAlgebra
    );

    for (int t = 0; t < 1000; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}




TEST(CPUMatrixTest, FulltestMatrixCPU) {

    Grid3D current = runCPU_Matrix();

    EXPECT_NEAR(current(3,4,5), 101.18867361286134, 1e-2);
    EXPECT_NEAR(current(7,6,8), 109.55658367148229, 1e-2);
    EXPECT_NEAR(current(1,3,1), 92.75879343854106, 1e-2);

    EXPECT_NEAR(current(0,0,0), 98.45116695237250, 1e-2);
    EXPECT_NEAR(current(2,2,2), 95.05961247275025, 1e-2);
    EXPECT_NEAR(current(4,4,4), 98.60626975111873, 1e-2);
    EXPECT_NEAR(current(6,6,6), 103.59085975524636, 1e-2);
    EXPECT_NEAR(current(8,8,8), 104.74529697665720, 1e-2);

    EXPECT_NEAR(current(9,1,7), 100.95294482977657, 1e-2);
    EXPECT_NEAR(current(5,9,3), 99.16204961395044, 1e-2);
    EXPECT_NEAR(current(8,2,9), 111.30994932204399, 1e-2);
}

#ifdef ENABLE_CUDA
#include "../include/cudaHeaders/solverCUDA.cuh"
#include "../include/cudaHeaders/linearAlgebra.cuh"
#include "../include/cudaHeaders/boundaryConditions.cuh"

Grid3D runGPU_Matrix() {

    Grid3D current(10,10,10,0.1);
    current.fill(50.0);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CUDA_MATRIX);

    Grid3D next = current;

    LinearAlgebraCUDA linAlgebraCUDA(50);
    BoundaryConditions bc(types, values);
    SimulationGlobals globs;

    HeatSolverCUDAMatrix solver(
        current,
        current.nx(),
        current.ny(),
        current.nz(),
        1e-4,
        current.dx(),
        globs.dt,
        globs.k,
        bc,
        linAlgebraCUDA
    );

    for (int t = 0; t < 1000; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    return current;
}



TEST(CUDAMatrixTest, FulltestMatrixGPU) {

    Grid3D current = runGPU_Matrix();

    EXPECT_NEAR(current(3,4,5), 101.18867361286134, 1e-2);
    EXPECT_NEAR(current(7,6,8), 109.55658367148229, 1e-2);
    EXPECT_NEAR(current(1,3,1), 92.75879343854106, 1e-2);

    EXPECT_NEAR(current(0,0,0), 98.45116695237250, 1e-2);
    EXPECT_NEAR(current(2,2,2), 95.05961247275025, 1e-2);
    EXPECT_NEAR(current(4,4,4), 98.60626975111873, 1e-2);
    EXPECT_NEAR(current(6,6,6), 103.59085975524636, 1e-2);
    EXPECT_NEAR(current(8,8,8), 104.74529697665720, 1e-2);

    EXPECT_NEAR(current(9,1,7), 100.95294482977657, 1e-2);
    EXPECT_NEAR(current(5,9,3), 99.16204961395044, 1e-2);
    EXPECT_NEAR(current(8,2,9), 111.30994932204399, 1e-2);
}

TEST(CompareCPUGPUMatrixTest, CPUGPUEqualityMatrix) {

    Grid3D cpu = runCPU_Matrix();
    Grid3D gpu = runGPU_Matrix();

    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < 10; ++k)
                EXPECT_NEAR(cpu(i,j,k), gpu(i,j,k), 1e-2);
}

#endif