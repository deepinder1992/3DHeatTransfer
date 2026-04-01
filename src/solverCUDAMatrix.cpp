#include "solverCUDA.hpp"
#include <stdexcept>
#include <cuda_runtime.h>
#include "kernel.cuh"
#include "heatMatrixBuilder.hpp"


//Implict Matrix Solver
HeatSolverCUDAMatrix::HeatSolverCUDAMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k,
                                            const BoundaryConditions& bc, const LinearAlgebra& linAlgebra):
                                            A_(grid.totalCellsInGeometry()),alpha_(alpha), dx_(dx), dt_(dt), cond_(k),
                                            linAlgebra_(linAlgebra){
    assert(alpha> 0.0);
    assert(dx > 0.0);
    assert (dt > 0.0);
    coeff_ = alpha_*dt_/(dx_*dx_);

    A_ = implicitMatrix(grid, nx, ny, nz, coeff_, bc); //heat matrix builder
}


void HeatSolverCUDAMatrix::step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs,const BoundaryConditions& bc){

    size_type N = current.totalCellsInGeometry();
    assert(N == A_.rows());
    std::vector<double> b(N, 0.0);
    std::size_t counter = 0;
    
    for (auto& cell : current.activeIndices()) {
        auto [i,j,k] = cell;
        b[counter] = current(i,j,k);
        ++counter;
    }
    
    bc.applyBCsToRhsMatrix(current, current.nx(), current.ny(), current.nz(), dx_,
                            coeff_, cond_, b);
                             
    std::vector<double> x(N,0.0);
    

    
    linAlgebra_.conjugateGradientCUDA(A_, b, x, globs);

    counter = 0;
    for (auto& cell : current.activeIndices()) {
        auto [i,j,k] = cell;
        next(i,j,k) = x[counter];
        ++counter;
    }
}
