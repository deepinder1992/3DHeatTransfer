#include "linearAlgebra.hpp"
#include "kernel.cuh"


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



void conjugateGradientCUDA(const SparseMatrix& A, const std::vector<double>& b,
                            std::vector<double>& x, const SimulationGlobals& globs)
{
    int N = b.size();
    std::vector<double> r(N), p(N), Ap(N);

    std::size_t aValSize = A.values().size();
    std::size_t aRowPtrSize = A.rowPtr().size();
    std::size_t aColIdxSize = A.colIndex().size();
    dim3 blockDim1D(globs.blockDimX*globs.blockDimY*globs.blockDimZ);
    dim3 gridDims1D((N+globs.blockDimX*globs.blockDimY*globs.blockDimZ-1)
                        /(globs.blockDimX*globs.blockDimY*globs.blockDimZ));

    double* hostSum = new double;
    double* hostSum2 = new double;
    double alpha  = 1.0;


    ::allocateMemory(devSparseMatValues, devMemSpMatVals, aValSize);
    ::allocateMemory(devSparseMatRowPtr, devMemSpMatRowPtr, aRowPtrSize);
    ::allocateMemory(devSparseMatCols, devMemSpMatCols, aColIdxSize);
    ::allocateMemory(devBVector, devMemBVector, N);
    ::allocateMemory(devXVector, devMemXVector, N);
    ::allocateMemory(devApVector, devMemApVector, N);
    ::allocateMemory(devRVector, devMemRVector, N);
    ::allocateMemory(devBlockSums, devMemBlockSums, gridDims1D.x);
    ::allocateMemory(devSum, devMemSum, 1);
    ::allocateMemory(devSum2, devMemSum2, 1);


    cudaMemcpy(devSparseMatValues, A.values(), aValSize*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devSparseMatRowPtr, A.rowPtr(), aRowPtrSize*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devSparseMatCols, A.colIndex(), aColIdxSize*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devBVector, b.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devXVector, x.data(), N*sizeof(double), cudaMemcpyHostToDevice);



    sparseMultiply<<<gridDim1D, blockDim1D>>>(devSparseMatValues, devSparseMatCols, devSparseMatRowPtr, devXVector, devApVector, N);
    subtract<<<gridDim1D, blockDim1D>>>(devBVector, devApVector, devRVector, N);

    long int sharedMemSizeDot = blockDims1D.x;
    dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devRVector, devRVector, devBlockSums, N);

    long int sharedMemSizeSum = blockDims1D.x*sizeof(double);
    
    arrayAtomicAdd <<<gridDim1D, blockDim1D,sharedMemSizeSum>>>(devBlockSums, devSum, gridDims1D.x);
    cudaMemcpy(hostSum, devSum , sizeof(double), cudaMemcpyDeviceToHost);


    for (int iter = 0; iter < globs.maxIters; ++iter)
    {
        sparseMultiply<<<gridDim1D, blockDim1D>>>(devSparseMatValues, devSparseMatCols, devSparseMatRowPtr, devPVector, devApVector, N);
        //SparseMultiply(A, p.data(), Ap.data());
        long int sharedMemSizeDot = blockDims1D.x;
        dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devPVector, devAPVector, devBlockSums, N);
        long int sharedMemSizeSum = blockDims1D.x*sizeof(double);
        arrayAtomicAdd <<<gridDim1D, blockDim1D,sharedMemSizeSum>>>(devBlockSums, devSum2, gridDims1D.x);
        cudaMemcpy(hostSum2, devSum2 , sizeof(double), cudaMemcpyDeviceToHost);
    
        alpha = (*hostSum)/(*hostSum2);

        addSubtract<<<gridDim1D,blockDim1D>>>(devXVector , devPVector, devXVector alpha, N, 1.0);
        addSubtract<<<gridDim1D,blockDim1D>>>(devRVector , devAPVector, devRVector alpha, N, -1.0);

        long int sharedMemSizeDot = blockDims1D.x;
        dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devRVector, devRVector, devBlockSums, N);
        long int sharedMemSizeSum = blockDims1D.x*sizeof(double);
        arrayAtomicAdd <<<gridDim1D, blockDim1D,sharedMemSizeSum>>>(devBlockSums, devSum2, gridDims1D.x);
        cudaMemcpy(hostSum2, devSum2 , sizeof(double), cudaMemcpyDeviceToHost);

        double sqrtRsnew = std::sqrt(*hostSum2);
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<sqrtRsnew << std::endl;
            ++globs.totalIters;}
        if (sqrtRsnew < globs.tol) break;
        addSubtract<<<gridDim1D,blockDim1D>>>(devRVector , devPVector,devPVector, (*hostSum2/ (*hostSum)), N, 1.0);

        *hostSum = *hostSum2;
    }

}
// void launchdotBlock(double* a, double* b, double* blockSum, int N)
// {
//     dotBlock (a, b, blockSum, N);
// }
