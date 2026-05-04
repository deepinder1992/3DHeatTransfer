#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/voxelReader.hpp"
#include<iostream>


TEST(VoxelTest, UnableToFindFile) {
    Grid3D grid(50, 50, 50);
    EXPECT_THROW(VoxelReader("", grid), std::runtime_error);
};

TEST(VoxelTest, ReadCubeStl) {
    Grid3D grid(50, 50, 50);
    std::string stlFilePathStr = "../stlFiles/cube/cube.stl";
    VoxelReader(stlFilePathStr, grid);
    
    EXPECT_EQ(grid.boundaryIndices().size(), 14408);
    EXPECT_EQ(grid.interiorIndices().size(), 110592);
    EXPECT_EQ(grid.activeIndices().size(), 125000);

    EXPECT_EQ(grid.cellType(5,6,5), CellType::INTERIOR);
    EXPECT_EQ(grid.cellType(0,0,0), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(49,6,5), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(0,6,5), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(5,0,5), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(5,6,0), CellType::BOUNDARY);

    EXPECT_EQ(grid.faceType(0,10,15), FaceType::WALL);
    EXPECT_EQ(grid.faceType(15,0,15), FaceType::WALL);
    EXPECT_EQ(grid.faceType(15,15,0), FaceType::OUTLET);

    Vector result = grid.cellFaceNormalized(10, 0, 15);

    EXPECT_NEAR(result.x, 0.0, 1e-3);

    EXPECT_NEAR(result.y, -1.0, 1e-3);

    EXPECT_NEAR(result.z, 0.0, 1e-3);

    Vector result1 = grid.cellFaceNormalized(0, 10, 15);

    EXPECT_NEAR(result1.x, -1.0, 1e-3);

    EXPECT_NEAR(result1.y, 0.0, 1e-3);

    EXPECT_NEAR(result1.z, 0.0, 1e-3);
};

TEST(VoxelTest, ReadCyclinderStl) {
    Grid3D grid(50, 50, 50);
    std::string stlFilePathStr = "../stlFiles/cylinder/cylinder.stl";
    VoxelReader(stlFilePathStr, grid);
    
    EXPECT_EQ(grid.boundaryIndices().size(), 10640);
    EXPECT_EQ(grid.interiorIndices().size(), 87360);
    EXPECT_EQ(grid.activeIndices().size(), 98000);

    EXPECT_EQ(grid.cellType(25,25,25), CellType::INTERIOR);
    EXPECT_EQ(grid.cellType(49,25,25), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(0,25,25), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(0,22,20), CellType::BOUNDARY);

    EXPECT_EQ(grid.faceType(25,25,49), FaceType::INLET);
    EXPECT_EQ(grid.faceType(15,0,15), FaceType::NONE);
    EXPECT_EQ(grid.faceType(25,25,0), FaceType::OUTLET);

    Vector result = grid.cellFaceNormalized(25, 25, 0);

    EXPECT_NEAR(result.x, 0.0, 1e-3);
    EXPECT_NEAR(result.y, 0.0, 1e-3);
    EXPECT_NEAR(result.z, -1.0, 1e-3);

    Vector result1 = grid.cellFaceNormalized(49, 25, 25);

    EXPECT_NEAR(result1.x, 1.0, 1e-3);
    EXPECT_NEAR(result1.y, 0.04361913353204727, 1e-3);
    EXPECT_NEAR(result1.z, 0.0, 1e-3);
};


TEST(VoxelTest, ReadSemiCyclinderStl) {
    Grid3D grid(50, 50, 50);
    std::string stlFilePathStr =  "../stlFiles/semiCylinder/semicylinder.stl";
    VoxelReader(stlFilePathStr, grid);
    
    EXPECT_EQ(grid.boundaryIndices().size(), 7470);
    EXPECT_EQ(grid.interiorIndices().size(), 40080);
    EXPECT_EQ(grid.activeIndices().size(), 47550);

    EXPECT_EQ(grid.cellType(25,25,25), CellType::INTERIOR);
    EXPECT_EQ(grid.cellType(23,25,0), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(0,25,25), CellType::SOLID);
    EXPECT_EQ(grid.cellType(23,22,49), CellType::BOUNDARY);

    EXPECT_EQ(grid.faceType(25,25,49), FaceType::INLET);
    EXPECT_EQ(grid.faceType(15,0,15), FaceType::NONE);
    EXPECT_EQ(grid.faceType(25,25,0), FaceType::OUTLET);

    Vector result = grid.cellFaceNormalized(25, 25, 0);

    EXPECT_NEAR(result.x, 0.0, 1e-3);

    EXPECT_NEAR(result.y, 0.0, 1e-3);

    EXPECT_NEAR(result.z, -1.0, 1e-3);

    Vector result1 = grid.cellFaceNormalized(25, 25, 49);

    EXPECT_NEAR(result1.x, 0.0, 1e-3);

    EXPECT_NEAR(result1.y, 0.0, 1e-3);

    EXPECT_NEAR(result1.z, 1.0, 1e-3);
};


TEST(VoxelTest, ReadLChannelStl) {
    Grid3D grid(50, 50, 50);
    std::string stlFilePathStr =  "../stlFiles/L_Channel/l.stl";
    VoxelReader(stlFilePathStr, grid);
    
    EXPECT_EQ(grid.boundaryIndices().size(), 13110);
    EXPECT_EQ(grid.interiorIndices().size(), 80640);
    EXPECT_EQ(grid.activeIndices().size(), 93750);

    EXPECT_EQ(grid.cellType(25,25,25), CellType::SOLID);
    EXPECT_EQ(grid.cellType(23,25,0), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(0,25,25), CellType::BOUNDARY);
    EXPECT_EQ(grid.cellType(3,3,45), CellType::INTERIOR);

    EXPECT_EQ(grid.faceType(5,0,7), FaceType::WALL);
    EXPECT_EQ(grid.faceType(15,0,15), FaceType::WALL);
    EXPECT_EQ(grid.faceType(10,49,10), FaceType::OUTLET);

    Vector result = grid.cellFaceNormalized(5, 5, 0);

    EXPECT_NEAR(result.x, 0.0, 1e-3);

    EXPECT_NEAR(result.y, 0.0, 1e-3);

    EXPECT_NEAR(result.z, -1.0, 1e-3);

    Vector result1 = grid.cellFaceNormalized(23, 25, 0);

    EXPECT_NEAR(result1.x, 0.0, 1e-3);

    EXPECT_NEAR(result1.y, 0.0, 1e-3);

    EXPECT_NEAR(result1.z, -1.0, 1e-3);
};