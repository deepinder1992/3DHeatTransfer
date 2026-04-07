#pragma once
#include <cmath>
#include "sparseMatrix.hpp"
#include <iostream>
#include "simGlobals.hpp"
#include <cuda_runtime.h>


class LinearAlgebra{
    public:         
         LinearAlgebra(int maxIters):_maxIters(maxIters){};
        ~LinearAlgebra() { if (devSparseMatValues) cudaFree(devSparseMatValues);
                           if (devSparseMatRowPtr) cudaFree(devSparseMatRowPtr);
                           if (devSparseMatCols)   cudaFree(devSparseMatCols);
                           if (devXVector) cudaFree(devXVector);
                           if (devBVector) cudaFree(devBVector);
                           if (devRVector) cudaFree(devRVector);
                           if (devPVector) cudaFree(devPVector);
                           if (devSum) cudaFree(devSum);
                           if (devSum2) cudaFree(devSum2);
                           if (devBlockSums) cudaFree(devBlockSums);}

        double dot (const std::vector<double>& a, const std::vector<double>& b);

        void SparseMultiply (const SparseMatrix& A, const double* x, double* y);

        void conjugateGradient(const SparseMatrix& A, const std::vector<double>& b,
                                std::vector<double>& x, const SimulationGlobals& globs);

        void implicitJacobiCPU(size_type nx, size_type ny, size_type nz, const double coeff_, double& maxerror,
                                 Grid3D* oldGrid, Grid3D* newGrid, const Grid3D& current);
        //cuda functions
        void implicitJacobiCUDA(double* oldVal, double* newVal, double* currentVal,  std::size_t (*intIndices)[3], std::size_t nIntIdxs,
                                         std::size_t nx, std::size_t ny, std::size_t nz, double coeff_, dim3 grid, dim3 block);

        void maxErrorCUDA(double* oldVal, double* newVal,  double* maxBlockError, std::size_t (*intIndices)[3], std::size_t nIntIdxs, std::size_t nx, std::size_t ny,
                    dim3 grid, dim3 block, size_t sharedMemSize);
        
        void conjugateGradientCUDA(const SparseMatrix& A, const std::vector<double>& b,
                            std::vector<double>& x, const SimulationGlobals& globs);

        int maxIters() const noexcept{return _maxIters;}

        //in case of non convergnce increase max iters but cap at 2000 to avoid infinte loop
        void adjustMaxItersIfNeeded(int iter){
                    if (iter==(_maxIters-1) && _maxIters<2000){
                        _maxIters = static_cast<int> (_maxIters*1.5);
                        std::cout<< "Inner Solution did not converge increasing maxIters to " << _maxIters<<std::endl;
                    } }
                        


private:
    double *devSparseMatValues = nullptr, *devBVector = nullptr, *devXVector = nullptr,
           *devApVector = nullptr, *devRVector = nullptr, *devPVector = nullptr,
           *devSum = nullptr, *devSum2 = nullptr, *devBlockSums = nullptr;

    std::size_t *devSparseMatRowPtr = nullptr, *devSparseMatCols = nullptr;

    size_type devMemSpMatVals = 0, devMemSpMatRowPtr = 0, devMemSpMatCols = 0,
              devMemBVector = 0, devMemXVector = 0, devMemApVector = 0,
              devMemRVector = 0, devMemBlockSums = 0, devMemSum = 0,
              devMemSum2 = 0, devMemPVector = 0;

    int _maxIters = 0;
};