#include "solverCUDA.hpp"
#include <stdexcept>
#include <cuda_runtime.h>
#include "kernel.cuh"
#include "heatMatrixBuilder.hpp"


//Implict Matrix Solver
// HeatSolverCUDAMatrix::HeatSolverCUDAMatrix(size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k,
//                                             const BoundaryConditions& bc, const LinearAlgebra& linAlgebra):
//                                             A_(nx*ny*nz), alpha_(alpha), dx_(dx), dt_(dt), cond_(k),
//                                             linAlgebra_(linAlgebra){
//     assert(alpha> 0.0);
//     assert(dx > 0.0);
//     assert (dt > 0.0);
//     coeff_ = alpha_*dt_/(dx_*dx_);

//     //A_ = implicitMatrix(nx, ny, nz, coeff_, bc); //heat matrix builder
// }


void HeatSolverCUDAMatrix::step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs,const BoundaryConditions& bc){
    assert(current.nx() == next.nx());
    assert(current.ny() == next.ny());
    assert(current.nz() == next.nz());
    size_type N = current.size();

    std::vector<double> b(current.data(), current.data()+N);
    // bc.applyBCsToRhsMatrix(current.nx(), current.ny(), current.nz(),
    //                             dx_, coeff_, cond_, b);
                             
    std::vector<double> x(N,0.0);
    
    linAlgebra_.conjugateGradientCUDA(A_, b, x, globs);

    for (size_type i = 0; i < N ; ++i){
        next.data()[i] = x[i];
    }

}
