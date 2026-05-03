#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include<iostream>

TEST(GridTest, BasicTest) {
    Grid3D g(10, 12, 8, 0.1);
    EXPECT_EQ(g.nx(), 10);
    EXPECT_EQ(g.ny(), 12);
    EXPECT_EQ(g.nz(), 8);
    EXPECT_EQ(g.size(), 960);
};

TEST(GridTest, FillFunction) {
    Grid3D g(10, 12, 8, 0.1);
    g.fill(42.5);
    EXPECT_EQ(g(5,6,3), 42.5);

};

TEST(GridTest, IndexFuncs) {
    Grid3D g(10, 12, 8, 0.1);
    g(5,6,5) = 44;
    EXPECT_EQ(g(5,6,5), 44);

    g.cellType(5,6,5) = CellType::BOUNDARY;
    EXPECT_EQ(g.cellType(5,6,5), CellType::BOUNDARY);

    g.faceType(5,6,5) = FaceType::INLET;
    EXPECT_NE(g.faceType(5,6,5), FaceType::WALL);

    g.cellFaceNormal(1,3,5) = {1.0, 1.0, 1.0}; 

    Vector result = g.cellFaceNormalized(1, 3, 5);

    EXPECT_NEAR(result.x, 0.577, 1e-3) 
        << "Expected x = 0.577, but got x = " << result.x;

    EXPECT_NEAR(result.y, 0.577, 1e-3) 
        << "Expected y = 0.577, but got y = " << result.y;

    EXPECT_NEAR(result.z, 0.577, 1e-3) 
        << "Expected z = 0.577, but got z = " << result.z;
};

TEST(GridTest, Constructor_ThrowsOnInvalidDimensions) {
    EXPECT_THROW(Grid3D g(0, 10, 10, 0.1), std::invalid_argument);
    EXPECT_THROW(Grid3D g(10, 10, 10, 0.0), std::invalid_argument);
    EXPECT_THROW(Grid3D g(-10, 10, 10, 0.0), std::length_error);
    EXPECT_THROW(Grid3D g(10, 10, 10, -0.1), std::invalid_argument);
}
