#pragma once
#include <cmath>
#include "sparseMatrix.hpp"
#include <iostream>
#include "simGlobals.hpp"
#include <cuda_runtime.h>


class LinearAlgebra{
    public:         
        //LinearAlgerbra();
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

        //cuda functions

        void implicitJacobiCPU(size_type nx, size_type ny, size_type nz, const double coeff_, double& maxerror,
                                 Grid3D* oldGrid, Grid3D* newGrid, const Grid3D& current);
            
        void implicitJacobiCUDA(double* oldVal, double* newVal, double* currentVal, int nx, int ny, int nz,
                                    double coeff_, dim3 grid, dim3 block);

        void maxErrorCUDA(double* oldVal, double* newVal,  double* maxBlockError, int N, int nx, int ny,
                            dim3 grid, dim3 block, size_t sharedMemSize);
        
        void conjugateGradientCUDA(const SparseMatrix& A, const std::vector<double>& b,
                            std::vector<double>& x, const SimulationGlobals& globs);

    private:
        double* devSparseMatValues = nullptr;
        std::size_t* devSparseMatRowPtr = nullptr;
        std::size_t* devSparseMatCols   = nullptr;
        double* devBVector = nullptr;
        double* devXVector = nullptr;
        double* devApVector = nullptr;
        double* devRVector = nullptr;
        double* devPVector = nullptr;
        double* devSum = nullptr;
        double* devSum2 = nullptr;
        double* devBlockSums = nullptr;

        size_type devMemSpMatVals = 0;
        size_type devMemSpMatRowPtr = 0;
        size_type devMemSpMatCols = 0;
        size_type devMemBVector = 0;
        size_type devMemXVector = 0;
        size_type devMemApVector = 0;
        size_type devMemRVector = 0;
        size_type devMemBlockSums = 0;
        size_type devMemSum = 0;
        size_type devMemSum2 = 0;
        size_type devMemPVector = 0;

};