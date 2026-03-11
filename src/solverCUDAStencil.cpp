#include "solverCUDA.hpp"
#include <stdexcept>
#include <cuda_runtime.h>
#include "kernel.cuh"

HeatSolverCUDAStencil::HeatSolverCUDAStencil(double alpha, double dx, double dt, const LinearAlgebra& linAlgebra):
                                alpha_(alpha), dx_(dx),dt_(dt),linAlgebra_(linAlgebra)
                            {  assert(alpha> 0.0);
                               assert(dx > 0.0);
                               assert (dt > 0.0);
                               coeff_ = alpha_*dt_/(dx_*dx_);
                        };

void HeatSolverCUDAStencil::step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc){

    size_type N = current.size();
    const std::size_t nx = current.nx();
    const std::size_t ny = current.ny();
    const std::size_t nz = current.nz();

    ::allocateMemory(devCurrent, devMemCurrGrdSize, N);
    ::allocateMemory(devNext, devMemNextGrdSize, N);
    ::allocateMemory(devOld, devMemOldGrdSize, N);
   
    cudaMemcpy(devCurrent, current.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devOld, next.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devNext, next.data(), N*sizeof(double), cudaMemcpyHostToDevice);

    
    dim3 blockDims(globs.blockDimX,globs.blockDimY,globs.blockDimZ);
    dim3 gridDims((nx+globs.blockDimX-1)/globs.blockDimX,
                    (ny+globs.blockDimY-1)/globs.blockDimY,
                    (nz+globs.blockDimZ-1)/globs.blockDimZ);


    size_type numBlocks = gridDims.x*gridDims.y*gridDims.z;

    ::allocateMemory(devMaxBlockError, devMemBlockErrorSize, numBlocks);
    
    for (int iter = 0; iter<linAlgebra_.maxIters();++iter){
        bc.applyBCsToStencilCUDA(devOld, dx_,nx, ny, nz, globs.k, gridDims, blockDims);   
        // cudaMemcpy(next.data(),devOld,N*sizeof(double),cudaMemcpyDeviceToHost);
        // bc.applyBCsToStencil(next, globs.dx, globs.k);
        // cudaMemcpy(devOld,next.data(),N*sizeof(double),cudaMemcpyHostToDevice);
        linAlgebra_.implicitJacobiCUDA(devOld, devNext, devCurrent, nx, ny, nz, coeff_, gridDims, blockDims);
        
        if (globs.verbosity & SimulationGlobals::VERB_MEDIUM){
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "CUDA error Jacobi: " << cudaGetErrorString(err) << std::endl;
            }
        }
        std::size_t sharedMemSize = blockDims.x*blockDims.y*blockDims.z*sizeof(double);
        linAlgebra_.maxErrorCUDA(devOld, devNext, devMaxBlockError, N, nx, ny, gridDims, blockDims, sharedMemSize);
        if (globs.verbosity & SimulationGlobals::VERB_MEDIUM){
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "CUDA error maxError Call: " << cudaGetErrorString(err) << std::endl;
            }
        }

        double hostMaxBlockError[numBlocks];
        cudaMemcpy(hostMaxBlockError, devMaxBlockError, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

        double maxErr = 0.0;
        for (std::size_t b = 0; b < numBlocks; ++b)
            maxErr = std::max(maxErr, hostMaxBlockError[b]);
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<maxErr<< std::endl;
            ++globs.totalIters;    
        }     
           
        if (maxErr<globs.tol)break;

        std::swap(devOld,devNext);
        linAlgebra_.adjustMaxItersIfNeeded(iter);
    }   
    cudaMemcpy(next.data(), devOld, N*sizeof(double), cudaMemcpyDeviceToHost);

    //bc.applyBCsToStencil(next, globs.dx, globs.k);   

    
}


