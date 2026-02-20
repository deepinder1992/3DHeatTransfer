#include "grid.hpp"
#include "solverCPU.hpp"
#include <iostream>
#include "boundaryConditions.hpp"
#include "outputWriter.hpp"
#include "simGlobals.hpp"


int main (){

    std::size_t nx = 31;
    std::size_t ny = 31;
    std::size_t nz = 31;

     
    SimulationGlobals globs;

    Grid3D current(nx, ny, nz);
    Grid3D next(nx, ny, nz);

    current.fill(75.0);

    
    BoundaryConditions bc(globs.types, globs.values);
    bc.applyBCsToStencil(current, globs.dx, globs.k);

    HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt);
    //HeatSolverCPUMatrix solver(nx, ny, nz, globs.alpha, globs.dx, globs.dt, globs.k, bc);

    BinaryWriter binWriter("../BinaryOutput", "temperature");
    VTKWriter vtkWriter("../VTKOutput", "temperature");
                                

    for (globs.t = 0; globs.t<globs.steps; ++globs.t){
        solver.step(current, next, globs, bc);        
        if (globs.t%globs.writeInterval == 0){
            binWriter.write(next, globs.t);
            vtkWriter.write(next, globs.t);
        }
        std::swap(current,next);

        std::cout << "Step "<< globs.t+1 <<": center "<< current(nx/2,ny/2,nz/2)<< std::endl;
    }

    std::cout << "Simulation Complete!" << std::endl;

    return 0;
}