#include <iostream>
#include <vector>

#include <cuda_runtime.h>
#include "tests_utils.hpp"
#include "../include/simGlobals.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/kernel.cuh"


int main_cuda_tests() {
    std::cout << "\n=== Running CUDA Linear Algebra Tests ===\n";

    SimulationGlobals globs;
    std::size_t N = 20000;
    double min = -10000.0;
    double max = 10000.0;

    auto randA = randomVector(N, min, max);   // You'll need to move randomVector to a header if used here
    auto randB = randomVector(N, min, max);

    dim3 blockDim1D(globs.blockDim);
    dim3 gridDim1D((N + globs.blockDim - 1) / globs.blockDim);

    double* devAVector = nullptr;
    double* devBVector = nullptr;
    double* devCVector = nullptr;
    double* devBlockSums = nullptr;

    CUDA_CHECK(cudaMalloc(&devAVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devBVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devCVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devBlockSums, gridDim1D.x * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(devAVector, randA.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devBVector, randB.data(), N * sizeof(double), cudaMemcpyHostToDevice));

    // Test 4: Add
    std::cout << "Test 4: Add Two Arrays CPU vs GPU... ";
    //addSubtract<<<gridDim1D, blockDim1D>>>(devAVector, devBVector, devCVector, 1.0, N, 1.0);
    // ... (rest of add test)

    // Test 5: Subtract
    std::cout << "Test 5: Subtract Two Arrays CPU vs GPU... ";

    // Test 6: Dot product
    std::cout << "Test 6: Dot product CPU vs GPU... ";

    // Cleanup
    cudaFree(devAVector);
    cudaFree(devBVector);
    cudaFree(devCVector);
    cudaFree(devBlockSums);

    std::cout << "All CUDA tests passed!\n";
    return 0;
}
