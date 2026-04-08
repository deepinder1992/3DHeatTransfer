#pragma once
#include <cuda_runtime.h>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <array>
#include "simGlobals.hpp"

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



__global__  void implicitJacobiKernel(double* oldVal, double* newVal, double* currentVal, std::size_t (*intIndices)[3], std::size_t nIntIdxs,
                                        std::size_t nx, std::size_t ny, std::size_t nz, double coeff_);

__global__ void addSubtract(double*a , double* b , double* c , double alpha, std::size_t N, double sign);

__global__ void dotBlock (const double* a, const double* b, double* blockSum, std::size_t N);

__global__ void sparseMultiply (const double* values, const std::size_t* cols, const std::size_t* rowPtr, const double*x, double* y, std::size_t N);

__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError,std::size_t (*intIndices)[3], std::size_t nIntIdxs, std::size_t nx, std::size_t ny);

__global__ void arraySumReduction (double* a, double* blockSum, std::size_t n);

__global__ void applyBCsToStencilKern(double* grid, double *oldGrid, std::size_t nx, std::size_t ny, std::size_t nz, double dx, 
                                    std::size_t (*bcIndices)[3], FaceType* faceTypes,std::size_t nBcCells,
                                    NeighbourType* nbrType, std::size_t* nbrOffset, float (*devCellNormals)[3], double  cond, const BCType types_[3],
                                    const double values_[3]);

double arraySum(double* d_a, std::size_t n);






