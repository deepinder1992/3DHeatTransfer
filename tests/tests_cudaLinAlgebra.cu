// tests/test_cuda_linearalgebra.cpp
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include "./include/tests_utils.hpp"
#include "../include/cudaHeaders/linearAlgebra.cuh"
#include "../include/cudaHeaders/kernel.cuh"
#include "../include/linearAlgebra.hpp"




TEST(CudaLinearAlgebraTest, AddTwoArrays_CPU_vs_GPU) {
    SimulationGlobals globs;
    const std::size_t N = 20000;
    const double min_val = -10000.0;
    const double max_val = 10000.0;

    auto randA = randomVector(N, min_val, max_val);
    auto randB = randomVector(N, min_val, max_val);

    dim3 blockDim1D(globs.blockDim);
    dim3 gridDim1D((N + globs.blockDim - 1) / globs.blockDim);

    double* devAVector = nullptr;
    double* devBVector = nullptr;
    double* devCVector = nullptr;

    CUDA_CHECK(cudaMalloc(&devAVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devBVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devCVector, N * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(devAVector, randA.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devBVector, randB.data(), N * sizeof(double), cudaMemcpyHostToDevice));

    // GPU Addition
    addSubtract<<<gridDim1D, blockDim1D>>>(devAVector, devBVector, devCVector, 1.0, N, 1.0);
    CUDA_CHECK(cudaGetLastError());

    std::vector<double> cudaRes(N);
    CUDA_CHECK(cudaMemcpy(cudaRes.data(), devCVector, N * sizeof(double), cudaMemcpyDeviceToHost));

    // CPU Addition
    std::vector<double> cpuRes(N);
    std::transform(randA.begin(), randA.end(), randB.begin(), cpuRes.begin(),
                   std::plus<double>());

    // Compare
    EXPECT_TRUE(std::equal(cpuRes.begin(), cpuRes.end(), cudaRes.begin(),
                [](double a, double b) { return std::abs(a - b) < 1e-10; }))
        << "GPU Add does not match CPU Add";

    // Cleanup
    CUDA_CHECK(cudaFree(devAVector));
    CUDA_CHECK(cudaFree(devBVector));
    CUDA_CHECK(cudaFree(devCVector));
}

TEST(CudaLinearAlgebraTest, SubtractTwoArrays_CPU_vs_GPU) {
    SimulationGlobals globs;
    const std::size_t N = 20000;

    auto randA = randomVector(N, -10000.0, 10000.0);
    auto randB = randomVector(N, -10000.0, 10000.0);

    dim3 blockDim1D(globs.blockDim);
    dim3 gridDim1D((N + globs.blockDim - 1) / globs.blockDim);

    double* devAVector = nullptr;
    double* devBVector = nullptr;
    double* devCVector = nullptr;

    CUDA_CHECK(cudaMalloc(&devAVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devBVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devCVector, N * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(devAVector, randA.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devBVector, randB.data(), N * sizeof(double), cudaMemcpyHostToDevice));

    // GPU Subtraction
    addSubtract<<<gridDim1D, blockDim1D>>>(devAVector, devBVector, devCVector, 1.0, N, -1.0);
    CUDA_CHECK(cudaGetLastError());

    std::vector<double> cudaRes(N);
    CUDA_CHECK(cudaMemcpy(cudaRes.data(), devCVector, N * sizeof(double), cudaMemcpyDeviceToHost));

    // CPU Subtraction
    std::vector<double> cpuRes(N);
    std::transform(randA.begin(), randA.end(), randB.begin(), cpuRes.begin(),
                   std::minus<double>());

    EXPECT_TRUE(std::equal(cpuRes.begin(), cpuRes.end(), cudaRes.begin(),
                [](double a, double b) { return std::abs(a - b) < 1e-10; }))
        << "GPU Subtract does not match CPU Subtract";

    CUDA_CHECK(cudaFree(devAVector));
    CUDA_CHECK(cudaFree(devBVector));
    CUDA_CHECK(cudaFree(devCVector));
}

TEST(CudaLinearAlgebraTest, DotProduct_CPU_vs_GPU) {
    SimulationGlobals globs;
    const std::size_t N = 20000;

    auto randA = randomVector(N, -10000.0, 10000.0);
    auto randB = randomVector(N, -10000.0, 10000.0);

    dim3 blockDim1D(globs.blockDim);
    dim3 gridDim1D((N + globs.blockDim - 1) / globs.blockDim);

    double* devAVector = nullptr;
    double* devBVector = nullptr;
    double* devBlockSums = nullptr;

    CUDA_CHECK(cudaMalloc(&devAVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devBVector, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&devBlockSums, gridDim1D.x * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(devAVector, randA.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(devBVector, randB.data(), N * sizeof(double), cudaMemcpyHostToDevice));

    // GPU Dot Product
    std::size_t sharedMemSize = blockDim1D.x * sizeof(double);
    dotBlock<<<gridDim1D, blockDim1D, sharedMemSize>>>(devAVector, devBVector, devBlockSums, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    double cudaDot = arraySum(devBlockSums, gridDim1D.x);

    // CPU Dot Product
    LinearAlgebra linAlg(globs.maxIters);
    double cpuDot = linAlg.dot(randA, randB);

    EXPECT_NEAR(cpuDot, cudaDot, 0.1)
        << "GPU Dot Product does not match CPU. CPU=" << cpuDot << ", GPU=" << cudaDot;

    // Cleanup
    CUDA_CHECK(cudaFree(devAVector));
    CUDA_CHECK(cudaFree(devBVector));
    CUDA_CHECK(cudaFree(devBlockSums));
}