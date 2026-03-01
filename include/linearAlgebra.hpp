#pragma once
#include <cmath>
#include "sparseMatrix.hpp"
#include <iostream>
#include "simGlobals.hpp"
#include <cuda_runtime.h>


class LinearAlgebra{
    public:         
        //LinearAlgerbra();
        //~LinearAlgebra();

        double dot (const std::vector<double>& a, const std::vector<double>& b);
        void SparseMultiply (const SparseMatrix& A, const double* x, double* y);
        void conjugateGradient(const SparseMatrix& A,
                       const std::vector<double>& b,
                       std::vector<double>& x,
                         const SimulationGlobals& globs);
        void implicitJacobiCPU(size_type nx, size_type ny, size_type nz, const double coeff_, double& maxerror,
                                 Grid3D* oldGrid, Grid3D* newGrid, const Grid3D& current);
            
        void implicitJacobiCUDA(double* oldVal, double* newVal, double* currentVal, int nx, int ny, int nz,
                          double coeff_, dim3 grid, dim3 block);

        void maxErrorCUDA(double* oldVal, double* newVal,  double* maxBlockError, int N, int nx, int ny,
                    dim3 grid, dim3 block, size_t sharedMemSize);
            
};