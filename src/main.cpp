#include "grid.hpp"
#include "solverCPU.hpp"
#include <iostream>
#include "boundaryConditions.hpp"

int main (){

    std::size_t nx = 10;
    std::size_t ny = 10;
    std::size_t nz = 10;


    double alpha = 0.01;
    double dx = 1.0;
    double dt = 0.1;

    Grid3D current(nx,ny,nz);
    Grid3D next(nx,ny,nz);

    current.fill(0.0);

    current(nx/2,ny/2,nz/2) = 100;

    HeatSolverCPU solver(alpha,dx,dt);

    std::array<BCType,6> types = {
                                    BCType::Dirichlet, BCType::Dirichlet, //xmin,max
                                    BCType::Dirichlet, BCType::Dirichlet,    //ymin,max
                                    BCType::Dirichlet, BCType::Dirichlet    //zmin,max
                                }

    std::array<BCType,6> values = {
                                    100,50,
                                    50,100,
                                    100,50
    
    
    BoundaryConditions bc(types,values);
                                }
    int steps = 5;

    for (int t = 0; t<steps; ++t){
        solver.step(current, next);
        bc.apply(next;)
        std::swap(current,next);

        std::cout << "Step "<<t+1<<": center "<< current(nx/2,ny/2,nx/2)<< std::endl;
    }

    std::cout << "Simulation Complete!" << std::endl;

    return 0;
}