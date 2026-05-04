#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/boundaryConditions.hpp"
#include "./include/tests_utils.hpp"
#include "../include/voxelReader.hpp"
#include<iostream>
#include <cmath>



void initSinusoidal(Grid3D& grid) {
    const double pi = 3.141592653589793;

    std::size_t nx = grid.nx();
    std::size_t ny = grid.ny();
    std::size_t nz = grid.nz();
    double dx = grid.dx();

    for (std::size_t k = 0; k < nz; ++k) {
        double z = k * dx;

        for (std::size_t j = 0; j < ny; ++j) {
            double y = j * dx;

            for (std::size_t i = 0; i < nx; ++i) {
                double x = i * dx;

                grid(i,j,k) =
                    std::sin(pi * x) *
                    std::sin(pi * y) *
                    std::sin(pi * z);
            }
        }
    }
}

double analytical(double x, double y, double z,
                  double t, double alpha) {
    const double pi = 3.141592653589793;

    return std::sin(pi * x) *
           std::sin(pi * y) *
           std::sin(pi * z) *
           std::exp(-3.0 * pi * pi * alpha * t);
}


TEST(CPUTestSinuSoidal, CPuMatrixSinuSoidal){
    //nx*nxy*nz 10*10*10 dx = 1/(nx-1)
    std::size_t nx = 50, ny = 50, nz = 50;

    Grid3D current(nx,ny,nz);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CPU_MATRIX);
    
    Grid3D next = current;
    LinearAlgebra linAlgebra(50);

    SimulationGlobals globs;
    globs.dt = 1.0;
    globs.types =  {
                    BCType::Dirichlet, //inlet
                    BCType::Dirichlet,    //outlet
                    BCType::Dirichlet       //wall
                };  
    
    globs.values = { 0.0, //inlet
                    0.0, //outlet
                    0.0}; //wall
    BoundaryConditions bc(globs.types, globs.values);

    const double pi = 3.141592653589793;

    initSinusoidal(current);

    HeatSolverCPUMatrix solver(current, current.nx(), current.ny(), current.nz(),  globs.alpha , current.dx(),  globs.dt,  globs.k, bc, linAlgebra );

    for (int t = 0; t < 10; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    // analytical solution
    double center = nx*current.dx()/2;
    double T_analytical =   analytical(center, center, center,10*globs.dt, globs.alpha);

    EXPECT_NEAR(current(nx/2,ny/2,nz/2), T_analytical, 5e-2);

}


TEST(CPUTestSinuSoidal, CPuStencilSinuSoidal){
    //nx*nxy*nz 10*10*10 dx = 1/(nx-1)
    std::size_t nx = 50, ny = 50, nz = 50;

    Grid3D current(nx,ny,nz);

    VoxelReader("../stlFiles/cube/cube.stl", current);
    current.constructNeigbourMap(SolverType::CPU_MATRIX);
    
    Grid3D next = current;
    LinearAlgebra linAlgebra(50);

    SimulationGlobals globs;
    globs.dt = 1.0;
    globs.types =  {
                    BCType::Dirichlet, //inlet
                    BCType::Dirichlet,    //outlet
                    BCType::Dirichlet       //wall
                };  
    
    globs.values = { 0.0, //inlet
                    0.0, //outlet
                    0.0}; //wall
    BoundaryConditions bc(globs.types, globs.values);

    const double pi = 3.141592653589793;

    initSinusoidal(current);

    HeatSolverCPUStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebra);

    for (int t = 0; t < 10; ++t) {
        solver.step(current, next, globs, bc);
        std::swap(current, next);
    }

    // analytical solution
    double center = nx*current.dx()/2;
    double T_analytical =   analytical(center, center, center,10*globs.dt, globs.alpha);

    EXPECT_NEAR(current(nx/2,ny/2,nz/2), T_analytical, 5e-2);

}
