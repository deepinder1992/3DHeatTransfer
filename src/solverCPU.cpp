#include "solverCPU.hpp"
#include <cassert>
#include "linearAlgebra.hpp"
#include "heatMatrixBuilder.hpp"

HeatSolverCPUStencil::HeatSolverCPUStencil(double alpha, double dx, double dt)
: alpha_(alpha), dx_(dx), dt_(dt)   
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

    Grid3D bufferGrid(nx,ny,nz);
    //bufferGrid.fill(0.0);

    // old/ new here just mean the internal prev and next itertions of this time step in for loop below,
    // current and next signify timesteps
    Grid3D* oldGrid = &bufferGrid;
    Grid3D* newGrid = &next;

    *oldGrid = next;

    // std::cout<< "Next1"<<next.data()<< std::endl;
    //         for (size_type k = 1; k < nz-1; ++k){
    //         for (size_type j = 1; j < ny-1; ++j ){
    //             for (size_type i = 1; i < nx-1 ; ++i){
    //                 std::cout << next(i,j,k) << " ";
    //             }}}
    
    for (int iter = 0; iter<=globs.maxIters;++iter){
        double maxErr = 0.0;
        #pragma omp parallel for collapse(3) reduction(max:maxErr)
        for (size_type k = 1; k < nz-1; ++k){
            for (size_type j = 1; j < ny-1; ++j ){
                for (size_type i = 1; i < nx-1 ; ++i){
                    const double rhs  = current(i,j,k);

                    const double sum = (*oldGrid)(i+1,j,k)+ (*oldGrid)(i-1,j,k)
                                                + (*oldGrid)(i,j+1,k)+ (*oldGrid)(i,j-1,k)
                                                    + (*oldGrid)(i,j,k+1)+ (*oldGrid)(i,j,k-1);
                    
                    const double newVal = (rhs + coeff_*sum)/(1+6*coeff_);
                    
                    maxErr = std::max(maxErr, std::abs(newVal-(*oldGrid)(i,j,k)));
                    (*newGrid)(i,j,k) = newVal;
                }
            }        
        }
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            #pragma omp critical
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<maxErr<< std::endl;}
        if (maxErr<globs.tol)break;
        std::swap(oldGrid,newGrid);
    }
    if (newGrid !=&next) next = *newGrid;
    //     std::cout<< "Next2"<<next.data()<< std::endl;
    //         for (size_type k = 1; k < nz-1; ++k){
    //         for (size_type j = 1; j < ny-1; ++j ){
    //             for (size_type i = 1; i < nx-1 ; ++i){
    //                 std::cout << next(i,j,k) << " ";
    //             }}}
    
    bc.applyBCsToStencil(next, globs.dx, globs.k);
        // std::cout<< "Next3"<<next.data()<< std::endl;
        //     for (size_type k = 1; k < nz-1; ++k){
        //     for (size_type j = 1; j < ny-1; ++j ){
        //         for (size_type i = 1; i < nx-1 ; ++i){
        //             std::cout << next(i,j,k) << " ";
        //         }}}
    
}   


//Implict Matrix Solver
HeatSolverCPUMatrix::HeatSolverCPUMatrix(size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k,
     const BoundaryConditions& bc):A_(nx*ny*nz), alpha_(alpha), dx_(dx), dt_(dt), cond_(k){
    assert(alpha> 0.0);
    assert(dx > 0.0);
    assert (dt > 0.0);
    coeff_ = alpha_*dt_/(dx_*dx_);

    A_ = implicitMatrix(nx, ny, nz, coeff_, bc);
}


void HeatSolverCPUMatrix::step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs,const BoundaryConditions& bc){
    assert(current.nx() == next.nx());
    assert(current.ny() == next.ny());
    assert(current.nz() == next.nz());
    size_type N = current.size();

    std::vector<double> b(current.data(), current.data()+N);
    bc.applyBCsToRhsMatrix(current.nx(),
                             current.ny(),
                              current.nz(),
                                dx_,
                                 coeff_,
                                  cond_,
                                    b);
                             
    std::vector<double> x(N,0.0);
    
    conjugateGradient(A_, b, x, globs);

    for (size_type i = 0; i < N ; ++i){
        next.data()[i] = x[i];
    }

}
