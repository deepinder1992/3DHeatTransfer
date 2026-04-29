#include <gtest/gtest.h>
#include "tests_utils.hpp"
#include "../include/grid.hpp"

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

TEST(GridTest, Index) {
    Grid3D g(10, 12, 8, 0.1);
    g(5,6,5) = 44;
    EXPECT_EQ(g(5,6,5), 44);

};

