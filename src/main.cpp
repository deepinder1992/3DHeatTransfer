#include "grid.hpp"
#include "solverCPU.hpp"
#include "solverCUDA.hpp"
#include "boundaryConditions.hpp"
#include "outputWriter.hpp"
#include "simGlobals.hpp"
#include <chrono>
#include <iostream>


int main (){
    SimulationGlobals globs;
    
    std::size_t nx = globs.nx;
    std::size_t ny = globs.ny;
    std::size_t nz = globs.nz;   
    
    Grid3D current(nx, ny, nz);
    Grid3D next(nx, ny, nz);

    current.fill(75.0);

    BoundaryConditions bc(globs.types, globs.values);
    bc.applyBCsToStencil(current, globs.dx, globs.k);

  // HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt);
  //  HeatSolverCPUMatrix solver(nx, ny, nz, globs.alpha, globs.dx, globs.dt, globs.k, bc);
  HeatSolverCUDAStencil solver(globs.alpha, globs.dx, globs.dt);

    BinaryWriter binWriter("../BinaryOutput", "temperature");
    VTKWriter vtkWriter("../VTKOutput", "temperature");
   
    double maxDiff = 0.0;                      
    size_type gridSize = current.size();

    auto start = std::chrono::high_resolution_clock::now();
    for (globs.t = 0; globs.t<globs.steps; ++globs.t){
        solver.step(current, next, globs, bc);        
        if (globs.t%globs.writeInterval == 0){
            binWriter.write(next, globs.t);
            vtkWriter.write(next, globs.t);
        }

        std::swap(current,next);

        std::cout << "Step "<< globs.t+1 <<": center "<< current(nx/2,ny/2,nz/2)<< ", Max Diff: "<< maxDiff<< std::endl;
        
        maxDiff = 0.0;
        for (size_type i = 0; i < gridSize; ++i)
            {
                double diff = std::abs(current.data()[i] - next.data()[i]);
                if (diff > maxDiff)
                    maxDiff = diff;
            }
        
        if (maxDiff < globs.globalTol) break;
    }
    std::cout << "Simulation Completed in "<< globs.t<< " iterations!" << std::endl;
    
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed =   std::chrono::duration<double>(end - start).count();
    std::cout <<solver.name()<< " Total simulation time "<< elapsed << " seconds" << std::endl;
    
    if (globs.verbosity & SimulationGlobals::VERB_HIGH){
        std::cout <<solver.name()<< " Totat Internal Iters "<< globs.totalIters<< std::endl;}

    return 0;
}