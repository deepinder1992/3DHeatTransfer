#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include "../include/linearAlgebra.hpp"
#include "../include/grid.hpp"
#include "./include/tests_utils.hpp"



TEST(CPULinearAlgebraTest, DotProduct ) {
    std::vector<double> a(100);
    std::vector<double> b(100);

    for (int i = 0; i<100; ++i){
        a[i] = i;
        b[i] = 100-i; }
    
    LinearAlgebra linAlgebra = LinearAlgebra(500);

    EXPECT_EQ(linAlgebra.maxIters(), 500);

    EXPECT_DOUBLE_EQ(linAlgebra.dot(a,b),166650.0);

}


TEST(CPULinearAlgebraTest, SparseMultiply ) {
    //create a 10 by 10 matrix of -1 2 1 entriens in each row, 2,1 in first -1, 2 in last

    SparseMatrix mat = makeTestMatrix10();
    std::vector<double> a(10);
    std::vector<double> b(10);

    for (int i = 0; i<10; ++i)
        a[i] = i;
    
    LinearAlgebra linAlg = LinearAlgebra(10);

    linAlg.SparseMultiply(mat,a.data(),b.data());

    EXPECT_DOUBLE_EQ(b[0], 1);
    EXPECT_DOUBLE_EQ(b[1], 4);
    EXPECT_DOUBLE_EQ(b[2], 6);
    EXPECT_DOUBLE_EQ(b[3], 8);
    EXPECT_DOUBLE_EQ(b[4], 10);
    EXPECT_DOUBLE_EQ(b[5], 12);
    EXPECT_DOUBLE_EQ(b[6], 14);
    EXPECT_DOUBLE_EQ(b[7], 16);
    EXPECT_DOUBLE_EQ(b[8], 18);
    EXPECT_DOUBLE_EQ(b[9], 10);
}

TEST(CPULinearAlgebraTest, ConjugateGradientTest ) {
    //create a 10 by 10 matrix of -1 2 1 entriens in each row, 2,1 in first -1, 2 in last

    SparseMatrix mat = makeTestMatrix10();
    SimulationGlobals globs;
    std::vector<double> a(10);
    std::vector<double> b(10);

    for (int i = 0; i<10; ++i)
        a[i] = i;
    
    LinearAlgebra linAlg = LinearAlgebra(10);

    linAlg.conjugateGradient(mat,a,b, globs);

    EXPECT_NEAR(b[0], -90.3781, 1e-3);
    EXPECT_NEAR(b[1], 24.1504, 1e-3);
    EXPECT_NEAR(b[2], 151.813, 1e-3);
    EXPECT_NEAR(b[3], -7.25714, 1e-3);
    EXPECT_NEAR(b[4], 136.936, 1e-3);
    EXPECT_NEAR(b[5], 54.7415, 1e-3);
    EXPECT_NEAR(b[6], -659.074, 1e-3);
    EXPECT_NEAR(b[7], 347.92, 1e-3);
    EXPECT_NEAR(b[8], 819.769, 1e-3);
    EXPECT_NEAR(b[9], -393.509, 1e-3);
}

TEST(CPULinearAlgebraTest, ImplicitJacobiTest ) {
    //test on grid of 10*10*10

    auto oldGrid = std::make_unique<Grid3D>(10,10,10,0.1); 
    auto  newGrid = std::make_unique<Grid3D>(10,10,10,0.1); 
    
    Grid3D current(10,10,10,0.1);
    double maxErr = 1e-3;

    (*oldGrid).fill(50.0);
    current.fill(45.0);
    

    LinearAlgebra linAlg = LinearAlgebra(10);
    linAlg.implicitJacobiCPU(10,10,10, 26, maxErr, oldGrid.get(), newGrid.get(), current);

    EXPECT_NEAR((*newGrid)(3,4,5), 49.968, 1e-3);
    EXPECT_NEAR((*newGrid)(7,6,8), 49.968, 1e-3);
    EXPECT_NEAR((*newGrid)(1,3,1), 49.968, 1e-3);
    EXPECT_NEAR(maxErr, 0.031, 1e-3);
}


