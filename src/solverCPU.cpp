#include "solverCPU.hpp"
#include <cassert>
#include "heatMatrixBuilder.hpp"


HeatSolverCPUStencil::HeatSolverCPUStencil(double alpha, double dx, double dt, const LinearAlgebra& linAlgebra)
: alpha_(alpha), dx_(dx), dt_(dt), linAlgebra_(linAlgebra) 
{
    
    assert(alpha> 0.0);
    assert(dx > 0.0);
    assert (dt > 0.0);

    coeff_ = alpha_*dt_/(dx_*dx_);
}


void HeatSolverCPUStencil::step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc){

    assert(current.nx() == next.nx());
    assert(current.ny() == next.ny());
    assert(current.nz() == next.nz());

    const std::size_t nx = current.nx();
    const std::size_t ny = current.ny();
    const std::size_t nz = current.nz();
   // const double dx = current.dx();
    
    //auto deep copy
    Grid3D bufferGrid = current;//(nx,ny,nz,dx);
    //bufferGrid.fill(0.0);

    // old/ new here just mean the internal prev and next itertions of this time step in for loop below,
    // current and next signify timesteps
    Grid3D* oldGrid = &bufferGrid;
    Grid3D* newGrid = &next;

    *oldGrid = next;
    
    for (int iter = 0; iter<linAlgebra_.maxIters();++iter){
        bc.applyBCsToStencil(*oldGrid, oldGrid->dx(), globs.k);
        double maxErr = 0.0;
        linAlgebra_.implicitJacobiCPU(nx, ny, nz, coeff_, maxErr, oldGrid, newGrid, current);

        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            #pragma omp critical
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<maxErr<< std::endl;
            ++globs.totalIters;    }
        if (maxErr<globs.tol)break;
        std::swap(oldGrid,newGrid);
        linAlgebra_.adjustMaxItersIfNeeded(iter);
    }
    if (newGrid !=&next) next = *newGrid;
    
    //bc.applyBCsToStencil(next, globs.dx, globs.k);    
}   


//Implict Matrix Solver
HeatSolverCPUMatrix::HeatSolverCPUMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k,
                                            const BoundaryConditions& bc, const LinearAlgebra& linAlgebra):
                                            A_(grid.totalCellsInGeometry()), alpha_(alpha), dx_(dx), dt_(dt), cond_(k),
                                            linAlgebra_(linAlgebra){
    assert(alpha> 0.0);
    assert(dx > 0.0);
    assert (dt > 0.0);
    coeff_ = alpha_*dt_/(dx_*dx_);

    A_ = implicitMatrix(grid, nx, ny, nz, coeff_, bc); //heat matrix builder
}


void HeatSolverCPUMatrix::step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs,const BoundaryConditions& bc){

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

    linAlgebra_.conjugateGradient(A_, b, x, globs);

    counter = 0;
    for (auto& cell : current.activeIndices()) {
        auto [i,j,k] = cell;
        next(i,j,k) = x[counter];
        ++counter;
    }

}
