#include "linearAlgebra.hpp"
#include "kernel.cuh"

void LinearAlgebra::implicitJacobiCPU(size_type nx, size_type ny, size_type nz, const double coeff_, double& maxErr,
                                     Grid3D* oldGrid, Grid3D* newGrid, const Grid3D& current){
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
}

void LinearAlgebra::implicitJacobiCUDA(double* oldVal, double* newVal, double* currentVal, int nx, int ny, int nz,
                          double coeff_, dim3 grid, dim3 block)
    {         
        implicitJacobiKernel<<<grid, block>>>(oldVal, newVal, currentVal, nx, ny, nz, coeff_);
        cudaDeviceSynchronize();
    }

void LinearAlgebra::maxErrorCUDA(double* oldVal, double* newVal,  double* maxBlockError, int N, int nx, int ny,
                    dim3 grid, dim3 block, size_t sharedMemSize)
    {
        maxError<<<grid, block, sharedMemSize>>>( oldVal, newVal, maxBlockError, N, nx, ny);
    }
