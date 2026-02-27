#pragma once
#include <cuda_runtime.h>
#include <stdexcept>
#include <cstddef>
#include <cmath>


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


void launchImplicitJacobi(double* oldVal,
                          double* newVal,
                          double* currentVal,
                          int nx, int ny, int nz,
                          double coeff_,
                          dim3 grid,
                          dim3 block);

void launchMaxError(double* oldVal,
                    double* newVal,
                    double* maxBlockError,
                    int N,
                    int nx,
                    int ny,
                    dim3 grid,
                    dim3 block,
                    size_t sharedMemSize);