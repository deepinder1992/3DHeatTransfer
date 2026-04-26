#include <iostream>
#include <vector>
#include "tests_utils.hpp"
#include "../include/cudaHeaders/linearAlgebra.cuh"
#include "../include/cudaHeaders/kernel.cuh"
#include "../include/linearAlgebra.hpp"


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
    //GPU addtion 
    addSubtract<<<gridDim1D, blockDim1D>>>(devAVector, devBVector, devCVector, 1.0, N, 1.0);
    std::vector<double> cudaRes(N);
    CUDA_CHECK(cudaMemcpy(cudaRes.data(), devCVector, N*sizeof(double), cudaMemcpyDeviceToHost));

    // cpu Addition
    std::vector<double> cpuRes(N);
    std::transform(randA.begin(), randA.end(), randB.begin(), cpuRes.begin(), [](double a, double b){return a+b;});

    //match the results
    bool ok = std::equal(cpuRes.begin(), cpuRes.end(), cudaRes.begin(), [](double a, double b){ return std::abs(a - b) < 1e-10; });
    if (ok) {std::cout << "PASS\n"; } 
    else {std::cout << "FAIL\n"; 
        return 1;}

    // Test 5: Subtract
    std::cout << "Test 5: Subtract Two Arrays CPU and GPU match... ";
    //GPU Subtraction
    addSubtract<<<gridDim1D, blockDim1D>>>(devAVector, devBVector, devCVector, 1.0, N, -1.0);
    CUDA_CHECK(cudaMemcpy(cudaRes.data(), devCVector, N*sizeof(double), cudaMemcpyDeviceToHost));

    // cpu subtract
    std::transform(randA.begin(), randA.end(), randB.begin(), cpuRes.begin(), [](double a, double b){return a-b;});

    //match the results
    ok = std::equal(cpuRes.begin(), cpuRes.end(), cudaRes.begin(), [](double a, double b){ return std::abs(a - b) < 1e-10; });
    if (ok) {std::cout << "PASS\n"; } 
    else {std::cout << "FAIL\n"; 
        return 1;}

    std::cout << "Test 6: Dot product Two Arrays CPU and GPU match... ";
    
    //GPU Dot
    std::size_t sharedMemSizeDot = blockDim1D.x*sizeof(double);
    dotBlock<<<gridDim1D, blockDim1D, sharedMemSizeDot>>>(devAVector, devBVector, devBlockSums, N);
    CUDA_CHECK(cudaGetLastError());     
    CUDA_CHECK(cudaDeviceSynchronize()); 
    double cudaDot = arraySum(devBlockSums, gridDim1D.x);

    // cpu dot
    LinearAlgebra linAlg(globs.maxIters);
    double cpuDot = linAlg.dot(randA,randB);
    
    CUDA_CHECK(cudaFree(devBVector));
    CUDA_CHECK(cudaFree(devAVector));
    CUDA_CHECK(cudaFree(devCVector));
    CUDA_CHECK(cudaFree(devBlockSums));

    //match the results
    if (std::abs(cpuDot - cudaDot) > 0.1) {
        std::cout << "FAIL (computed=" << cudaDot << ", expected≈" << cpuDot << ", diff = "<< std::abs(cpuDot - cudaDot)<<")\n";
        return 1;
    }
    std::cout << "PASS \n";

    // Cleanup
    cudaFree(devAVector);
    cudaFree(devBVector);
    cudaFree(devCVector);
    cudaFree(devBlockSums);

    std::cout << "All CUDA tests passed!\n";
    return 0;
}
