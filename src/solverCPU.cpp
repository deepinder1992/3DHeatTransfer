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


void HeatSolverCPUStencil::step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs){

    assert(current.nx() == next.nx());
    assert(current.ny() == next.ny());
    assert(current.nz() == next.nz());

    const std::size_t nx = current.nx();
    const std::size_t ny = current.ny();
    const std::size_t nz = current.nz();


    for (int k = 1; k < nz-1; ++k){
        for (int j = 1; j < ny-1; ++j ){
            for (int i = 1; i < nx-1 ; ++i){
                const double center  = current(i,j,k);

                const double laplacian = current(i+1,j,k)+current(i-1,j,k)
                                            + current(i,j+1,k)+current(i,j-1,k)
                                                + current(i,j,k+1)+current(i,j,k-1)
                                                    -6*center;
                
                next(i,j,k) = center + coeff_*laplacian;
            }
        }        
    }

}


//Implict Matrix Solver
HeatSolverCPUMatrix::HeatSolverCPUMatrix(size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt):A_(nx*ny*nz){
    assert(alpha> 0.0);
    assert(dx > 0.0);
    assert (dt > 0.0);

    double coeff = alpha*dt/(dx*dx);

    A_ = implicitMatrix(nx ,ny ,nz ,coeff);
}


void HeatSolverCPUMatrix::step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs){

    size_type N = current.size();

    std::vector<double> b(current.data(),current.data()+N);
    std::vector<double> x(N,0.0);
    
    conjugateGradient(A_ ,b ,x , 50, 1e-6, globs);

    for (size_type i = 0; i < N ; ++i){
        next.data()[i] = x[i];
    }

}
