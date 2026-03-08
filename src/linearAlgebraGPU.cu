#include "linearAlgebra.hpp"
#include "kernel.cuh"


void LinearAlgebra::implicitJacobiCUDA(double* oldVal, double* newVal, double* currentVal, std::size_t nx, std::size_t ny, std::size_t nz,
                          double coeff_, dim3 grid, dim3 block)
    {         
        implicitJacobiKernel<<<grid, block>>>(oldVal, newVal, currentVal, nx, ny, nz, coeff_);
        cudaDeviceSynchronize();
    }

void LinearAlgebra::maxErrorCUDA(double* oldVal, double* newVal,  double* maxBlockError, std::size_t N, std::size_t nx, std::size_t ny,
                    dim3 grid, dim3 block, size_t sharedMemSize)
    {
        maxError<<<grid, block, sharedMemSize>>>( oldVal, newVal, maxBlockError, N, nx, ny);
    }



void LinearAlgebra::conjugateGradientCUDA(const SparseMatrix& A, const std::vector<double>& b,
                            std::vector<double>& x, const SimulationGlobals& globs)
{
    std::size_t N = b.size();
    std::vector<double> r(N), p(N), Ap(N);

    std::size_t aValSize = A.values().size();
    std::size_t aRowPtrSize = A.rowPtr().size();
    std::size_t aColIdxSize = A.colIndex().size();
    dim3 blockDim1D(globs.blockDimX*globs.blockDimY*globs.blockDimZ);
    dim3 gridDim1D((N+globs.blockDimX*globs.blockDimY*globs.blockDimZ-1)
                        /(globs.blockDimX*globs.blockDimY*globs.blockDimZ));


    double fac  = 1.0;
    std::size_t sz = 1;
    std::size_t sharedMemSizeDot = blockDim1D.x*sizeof(double);
    

    ::allocateMemory(devSparseMatValues, devMemSpMatVals, aValSize);
    ::allocateMemory(devSparseMatRowPtr, devMemSpMatRowPtr, aRowPtrSize);
    ::allocateMemory(devSparseMatCols, devMemSpMatCols, aColIdxSize);
    ::allocateMemory(devBVector, devMemBVector, N);
    ::allocateMemory(devXVector, devMemXVector, N);
    ::allocateMemory(devApVector, devMemApVector, N);
    ::allocateMemory(devRVector, devMemRVector, N);
    ::allocateMemory(devPVector, devMemPVector, N);
    ::allocateMemory(devBlockSums, devMemBlockSums, gridDim1D.x);
    ::allocateMemory(devSum, devMemSum, sz);
    ::allocateMemory(devSum2, devMemSum2, sz);


    CUDA_CHECK(cudaMemcpy(devSparseMatValues, A.values().data(), aValSize*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devSparseMatRowPtr, A.rowPtr().data(), aRowPtrSize*sizeof(std::size_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devSparseMatCols, A.colIndex().data(), aColIdxSize*sizeof(std::size_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devBVector, b.data(), N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devXVector, x.data(), N*sizeof(double), cudaMemcpyHostToDevice));



    sparseMultiply<<<gridDim1D, blockDim1D>>>(devSparseMatValues, devSparseMatCols, devSparseMatRowPtr, devXVector, devApVector, N);
    CUDA_CHECK(cudaGetLastError());     
    CUDA_CHECK(cudaDeviceSynchronize()); 

    addSubtract<<<gridDim1D, blockDim1D>>>(devBVector, devApVector, devRVector, 1.0, N, -1.0);
    CUDA_CHECK(cudaGetLastError());     
    CUDA_CHECK(cudaDeviceSynchronize()); 
    
    dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devRVector, devRVector, devBlockSums, N);
    CUDA_CHECK(cudaGetLastError());     
    CUDA_CHECK(cudaDeviceSynchronize()); 
    
    double rsold = arraySum(devBlockSums, gridDim1D.x);

    CUDA_CHECK(cudaMemcpy(devPVector, devRVector, N * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy( r.data(), devRVector, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy( p.data(), devRVector, N*sizeof(double), cudaMemcpyHostToDevice));

    for (int iter = 0; iter < _maxIters; ++iter)
    {
        sparseMultiply<<<gridDim1D, blockDim1D>>>(devSparseMatValues, devSparseMatCols, devSparseMatRowPtr, devPVector, devApVector, N);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devPVector, devApVector, devBlockSums, N);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 

        fac = rsold / arraySum(devBlockSums, gridDim1D.x);

        addSubtract<<<gridDim1D,blockDim1D>>>(devXVector , devPVector, devXVector, fac, N, 1.0);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        addSubtract<<<gridDim1D,blockDim1D>>>(devRVector , devApVector, devRVector, fac, N, -1.0);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        
        dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devRVector, devRVector, devBlockSums, N);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        double rsnew = arraySum(devBlockSums, gridDim1D.x);
        double sqrtRsnew = std::sqrt(rsnew);

        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<sqrtRsnew << std::endl;
            ++globs.totalIters;}
        if (sqrtRsnew < globs.tol) break;

        addSubtract<<<gridDim1D,blockDim1D>>>(devRVector , devPVector,devPVector, (rsnew/ rsold), N, 1.0);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 

        rsold = rsnew;
        adjustMaxItersIfNeeded(iter);
    }
    CUDA_CHECK(cudaMemcpy( x.data(), devXVector, N*sizeof(double), cudaMemcpyDeviceToHost));

}
