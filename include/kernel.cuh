#pragma once
#include <cuda_runtime.h>
#include <stdexcept>


template<typename T>
inline void allocateMemory(T*& ptr,
                           std::size_t& allocatedSize,
                           std::size_t N)
{
    if (N != allocatedSize) {
        if (ptr)
            cudaFree(ptr);

        if (cudaMalloc(&ptr, N * sizeof(T)) != cudaSuccess)
            throw std::runtime_error("cudaMalloc Failed");

        allocatedSize = N;
    }
}

__global__  void implicitJacobiKernel(double* oldVal, double* newVal, double* currentVal, int nx, int ny, int nz, double coeff_);


__global__ void addSubtract(double*a , double* b , double alpha, int N, double sign);


__global__ void dotBlock (double* a, double* b, double* blockSum, int N);


__global__ void sparseMultiply (const double* values, const std::size_t* cols, const std::size_t* rowPtr, const double*x, double* y, int N);


__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError, int N, int nx, int ny);


__global__ void arraySumReduction (double* a, double* blockSum, int n);


__global__ void arrayAtomicAdd (double* a, double* result, int n);






