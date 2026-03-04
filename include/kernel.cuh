#pragma once
#include <cuda_runtime.h>
#include <stdexcept>
#include <iostream>
#include <cstdlib>

#define CUDA_CHECK(call)                                                \
do {                                                                    \
    cudaError_t err = call;                                             \
    if (err != cudaSuccess) {                                           \
        std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__   \
                  << " -> " << cudaGetErrorString(err) << std::endl;    \
        std::exit(EXIT_FAILURE);                                        \
    }                                                                   \
} while (0)

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


__global__ void addSubtract(double*a , double* b , double* c , double alpha, int N, double sign);


__global__ void dotBlock (const double* a, const double* b, double* blockSum, int N);


__global__ void sparseMultiply (const double* values, const std::size_t* cols, const std::size_t* rowPtr, const double*x, double* y, int N);


__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError, int N, int nx, int ny);


__global__ void arraySumReduction (double* a, double* blockSum, int n);


__global__ void arrayAtomicAdd (double* a, double* result, int n);

double arraySum(double* d_a, int n);






