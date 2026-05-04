#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/boundaryConditions.hpp"
#include "../include/voxelReader.hpp"
#include "./include/tests_utils.hpp"
#include<iostream>



TEST(CPUBoundaryConditionsTest, BasicTest) {
    SimulationGlobals globs;
    BoundaryConditions bc = BoundaryConditions(globs.types, globs.values);
    
    std::array<BCType,3> types_ = bc.types();
    std::array<double,3> values_ = bc.values();
    
    EXPECT_EQ(types_[0], BCType::Neumann);
    EXPECT_EQ(types_[1], BCType::Neumann);
    EXPECT_EQ(types_[2], BCType::Dirichlet );

    EXPECT_EQ(values_[0], 500.0);
    EXPECT_EQ(values_[1], -500.0);
    EXPECT_EQ(values_[2], 100.0);

}

TEST(CPUBoundaryConditionsTest, ApplyBCToStencil) {

    SimulationGlobals globs;
    BoundaryConditions bc = BoundaryConditions(globs.types, globs.values);

    Grid3D grid(50, 50, 50);
    std::string stlFilePathStr = "../stlFiles/cube/cube.stl";
    VoxelReader(stlFilePathStr, grid);
    grid.constructNeigbourMap(SolverType::CPU_MATRIX);

    bc.applyBCsToStencil(grid, 0.01, 100.0);

    EXPECT_NEAR(grid(5,6,3), 0, 1e-3);
    EXPECT_NEAR(grid(0,6,3), 100, 1e-3);
    EXPECT_NEAR(grid(5,0,3), 100, 1e-3);
    EXPECT_NEAR(grid(5,6,0), -0.1, 1e-3);

    EXPECT_NEAR(grid(0,0,0), 0.0333333, 1e-3);
    EXPECT_NEAR(grid(49,0,0), 0.0333333, 1e-3);
    EXPECT_NEAR(grid(0,49,0), 0.0333333, 1e-3);
    EXPECT_NEAR(grid(0,49,49), 49.9667, 1e-3);
}


TEST(CPUBoundaryConditionsTest, ApplyBCToRHS) {

    SimulationGlobals globs;
    BoundaryConditions bc = BoundaryConditions(globs.types, globs.values);

    Grid3D grid(50, 50, 50);
    grid.fill(50.0);
    std::string stlFilePathStr = "../stlFiles/cube/cube.stl";
    VoxelReader(stlFilePathStr, grid);
    grid.constructNeigbourMap(SolverType::CPU_STENCIL);
    std::vector<double> b(125000);
    for (int i=0; i<125000; ++i){
        b[i] = i/1000;
    }
    bc.applyBCsToRhsMatrix(grid,grid.nx(), grid.ny(), 0.1, 385, 100, b);

    EXPECT_NEAR(b[0], 385, 1e-3);
    EXPECT_NEAR(b[10000], 154010, 1e-3);
    EXPECT_NEAR(b[30000], 154030, 1e-3);
    EXPECT_NEAR(b[40000], 154040, 1e-3);
    EXPECT_NEAR(b[60000], 154060, 1e-3);
}

