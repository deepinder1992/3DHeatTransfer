#include "grid.hpp"
#include "solverCPU.hpp"
#include <iostream>
#include "boundaryConditions.hpp"
#include "outputWriter.hpp"
#include "simGlobals.hpp"


int main (){

    std::size_t nx = 50;
    std::size_t ny = 50;
    std::size_t nz = 50;

     
    SimulationGlobals globs;
    globs.writeInterval = 1000;
    globs.steps = 10000;
    globs.verbosity = SimulationGlobals::VERB_HIGH;

    Grid3D current(nx,ny,nz);
    Grid3D next(nx,ny,nz);

    current.fill(75.0);

    current(nx/2,ny/2,nz/2) = 100;

    //HeatSolverCPUStencil solver(alpha,dx,dt);
    HeatSolverCPUMatrix solver(nx,ny,nz,globs.alpha,globs.dx,globs.dt);
   
    std::array<BCType,6> types = {
                                    BCType::Dirichlet, BCType::Dirichlet, //xmin,max
                                    BCType::Dirichlet, BCType::Dirichlet,    //ymin,max
                                    BCType::Dirichlet, BCType::Dirichlet    //zmin,max
                                };

    std::array<double,6> values = {
                                    100,50,
                                    50,100,
                                    100,50};
    
    
    BoundaryConditions bc(types,values);

    BinaryWriter binWriter("../BinaryOutput","temperature");
    VTKWriter vtkWriter("../VTKOutput","temperature");
                                


    for (globs.t = 0; globs.t<globs.steps; ++globs.t){
        solver.step(current, next,globs);
        bc.apply(next, globs.t);
        if (globs.t%globs.writeInterval == 0){
            binWriter.write(next,globs.t);
            vtkWriter.write(next,globs.t);
        }
        std::swap(current,next);

        std::cout << "Step "<<globs.t+1<<": center "<< current(nx/2,ny/2,nz/2)<< std::endl;
    }

    std::cout << "Simulation Complete!" << std::endl;

    return 0;
}