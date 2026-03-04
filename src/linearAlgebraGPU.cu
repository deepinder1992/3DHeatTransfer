#include "linearAlgebra.hpp"
#include "kernel.cuh"

#include <numeric> 
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



void LinearAlgebra::conjugateGradientCUDA(const SparseMatrix& A, const std::vector<double>& b,
                            std::vector<double>& x, const SimulationGlobals& globs)
{
    int N = b.size();
    std::vector<double> r(N), p(N), Ap(N);

    std::size_t aValSize = A.values().size();
    std::size_t aRowPtrSize = A.rowPtr().size();
    std::size_t aColIdxSize = A.colIndex().size();
    dim3 blockDim1D(globs.blockDimX*globs.blockDimY*globs.blockDimZ);
    dim3 gridDim1D((N+globs.blockDimX*globs.blockDimY*globs.blockDimZ-1)
                        /(globs.blockDimX*globs.blockDimY*globs.blockDimZ));

    double* hostSum = new double;
    double* hostSum2 = new double;
    double fac  = 1.0;
    std::size_t sz = 1;
    long int sharedMemSizeDot = blockDim1D.x*sizeof(double);
    

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
    
    //std::vector<double> hosstSum(gridDim1D.x);
    double sum = arraySum(devBlockSums, gridDim1D.x);
    //CUDA_CHECK(cudaMemcpy(hosstSum.data(), devBlockSums , gridDim1D.x*sizeof(double), cudaMemcpyDeviceToHost));
    //double sum = std::accumulate(hosstSum.begin(), hosstSum.end(), 0.0);
    *hostSum = sum;
   // arrayAtomicAdd <<<gridDim1D, blockDim1D,sharedMemSizeSum>>>(devBlockSums, devSum, gridDim1D.x);
    //cudaDeviceSynchronize();
    //cudaMemcpy(hostSum, devSum , sizeof(double), cudaMemcpyDeviceToHost);
    CUDA_CHECK(cudaMemcpy(devPVector, devRVector, N * sizeof(double), cudaMemcpyDeviceToDevice));

    double rsold = sum;
    CUDA_CHECK(cudaMemcpy( r.data(), devRVector, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy( p.data(), devRVector, N*sizeof(double), cudaMemcpyHostToDevice));
    for (int iter = 0; iter < globs.maxIters; ++iter)
    {
        sparseMultiply<<<gridDim1D, blockDim1D>>>(devSparseMatValues, devSparseMatCols, devSparseMatRowPtr, devPVector, devApVector, N);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devPVector, devApVector, devBlockSums, N);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        sum = arraySum(devBlockSums, gridDim1D.x);
        //std::vector<double> hosstSum2(gridDim1D.x);
        //CUDA_CHECK(cudaMemcpy(hosstSum2.data(), devBlockSums , gridDim1D.x*sizeof(double), cudaMemcpyDeviceToHost));
        //sum = std::accumulate(hosstSum2.begin(), hosstSum2.end(), 0.0);
        *hostSum2 = sum;

        double alpha = rsold / (*hostSum2);

        addSubtract<<<gridDim1D,blockDim1D>>>(devXVector , devPVector, devXVector, fac, N, 1.0);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        addSubtract<<<gridDim1D,blockDim1D>>>(devRVector , devApVector, devRVector, fac, N, -1.0);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 


        
        dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devRVector, devRVector, devBlockSums, N);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 
        sum = arraySum(devBlockSums, gridDim1D.x);
        // arrayAtomicAdd <<<gridDim1D, blockDim1D,sharedMemSizeSum>>>(devBlockSums, devSum2, gridDim1D.x);
        // cudaDeviceSynchronize();
        // cudaMemcpy(hostSum2, devSum2 , sizeof(double), cudaMemcpyDeviceToHost);
        //std::vector<double> hosstSum3(gridDim1D.x);
        //CUDA_CHECK(cudaMemcpy(hosstSum3.data(), devBlockSums , gridDim1D.x*sizeof(double), cudaMemcpyDeviceToHost));
        //sum = std::accumulate(hosstSum3.begin(), hosstSum3.end(), 0.0);
        *hostSum2 = sum;
        double sqrtRsnew = std::sqrt(*hostSum2);
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<sqrtRsnew << std::endl;
            ++globs.totalIters;}
        if (sqrtRsnew < globs.tol) break;

        addSubtract<<<gridDim1D,blockDim1D>>>(devRVector , devPVector,devPVector, (*hostSum2/ (*hostSum)), N, 1.0);
        CUDA_CHECK(cudaGetLastError());     
        CUDA_CHECK(cudaDeviceSynchronize()); 

        *hostSum = *hostSum2;

    }
    CUDA_CHECK(cudaMemcpy( x.data(), devXVector, N*sizeof(double), cudaMemcpyDeviceToHost));

}
// void launchdotBlock(double* a, double* b, double* blockSum, int N)
// {
//     dotBlock (a, b, blockSum, N);
// }
