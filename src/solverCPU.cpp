#include "solverCPU.hpp"
#include <cassert>
#include "linearAlgebra.hpp"
#include "heatMatrixBuilder.hpp"


/////////////////////////
//////////////////////////
//////////////////////
#include "extra.hpp"

HeatSolverCPUStencil::HeatSolverCPUStencil(double alpha, double dx, double dt)
: alpha_(alpha), dx_(dx), dt_(dt)   
{
    
    assert(alpha> 0.0);
    assert(dx > 0.0);
    assert (dt > 0.0);

    coeff_ = alpha_*dt_/(dx_*dx_);
}


void HeatSolverCPUStencil::step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs){

    assert(current.nx() == next.nx());
    assert(current.ny() == next.ny());
    assert(current.nz() == next.nz());

    const std::size_t nx = current.nx();
    const std::size_t ny = current.ny();
    const std::size_t nz = current.nz();

    for (int iter = 0; iter<=globs.maxIters;++iter){
        double maxErr = 0.0;
        for (int k = 1; k < nz-1; ++k){
            for (int j = 1; j < ny-1; ++j ){
                for (int i = 1; i < nx-1 ; ++i){
                    const double rhs  = current(i,j,k);

                    const double sum = next(i+1,j,k)+next(i-1,j,k)
                                                + next(i,j+1,k)+ next(i,j-1,k)
                                                    + next(i,j,k+1)+ next(i,j,k-1);
                    
                    const double newVal = (rhs + coeff_*sum)/(1+6*coeff_);
                    
                    maxErr = std::max(maxErr, std::abs(newVal-next(i,j,k)));
                    next(i,j,k) = newVal;
                }
            }        
        }
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<maxErr<< std::endl;}
        if (maxErr<globs.tol)break;
    }
}   


//Implict Matrix Solver
HeatSolverCPUMatrix::HeatSolverCPUMatrix(size_type nx, size_type ny, size_type nz, const SimulationGlobals& globs):A_(nx*ny*nz){
    assert(globs.alpha> 0.0);
    assert(globs.dx > 0.0);
    assert (globs.dt > 0.0);

    A_ = implicitMatrix(nx ,ny ,nz ,globs);
}


void HeatSolverCPUMatrix::step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs){
    assert(current.nx() == next.nx());
    assert(current.ny() == next.ny());
    assert(current.nz() == next.nz());
    size_type N = current.size();

    std::vector<double> b(current.data(),current.data()+N);
    applyBoundaryConditions(current.nx(),
                             current.ny(),
                              current.nz(),
                                globs,
                             b);
                             
    std::vector<double> x(N,0.0);
    
    conjugateGradient(A_ ,b ,x ,globs);

    for (size_type i = 0; i < N ; ++i){
        next.data()[i] = x[i];
    }

}
