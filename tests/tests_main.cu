#include <iostream>
#include <cmath>
#include <cassert>
#include <random>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cuda_runtime.h>
#include "../include/grid.hpp"
#include "../include/simGlobals.hpp"
#include "../include/boundaryConditions.hpp"
#include "../include/solverCPU.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/kernel.cuh"


//random Vector Generator 
template<typename T> 
std::vector<T> randomVector(std::size_t num, T min, T max){
        std::vector<T> randVec(num);
        std::mt19937 gen(std::random_device{}());
        auto dist = [=](){if constexpr(std::is_integral_v<T>){
            return std::uniform_int_distribution<T> (min,max);
        }
        else if constexpr (std::is_floating_point_v<T>)
        {
            return  std::uniform_real_distribution<T>(min,max);   
        }}();
        
        for(auto &x:randVec){
            x = dist(gen);
        }
        return randVec;
}


// Helper: 1D Gaussian analytical (infinite domain)
double analytical_1d(double x, double t, double alpha) {
    if (t <= 0.0) return (std::abs(x) < 1e-10) ? 100.0 : 0.0;
    return 100.0 * std::exp(-x*x / (4.0 * alpha * t)) / std::sqrt(4.0 * M_PI * alpha * t);
}

bool approx_equal(double a, double b, double tol = 1e-4) {
    return std::abs(a - b) < tol;
}

int main() {
    std::cout << "=== Running simple manual tests ===\n\n";

    // Test 1: Grid basics
    {
        std::cout << "Test 1: Grid indexing & fill... ";
        Grid3D g(10, 12, 8);
        if (g.nx() != 10 || g.ny() != 12 || g.nz() != 8 || g.size() != 960) {
            std::cout << "FAIL\n";
            return 1;
        }
        g.fill(42.5);
        if (!approx_equal(g(5,6,4), 42.5)) {
            std::cout << "FAIL\n";
            return 1;
        }
        std::cout << "PASS\n";
    }

    // Test 2: CPU stencil – diffusion from center
    {
        std::cout << "Test 2: Stencil diffusion (qualitative)... ";
        SimulationGlobals globs;
        globs.nx = globs.ny = globs.nz = 21;
        globs.dx = globs.lx / globs.nx;
        globs.dt = 0.005;
        globs.alpha = 0.1;
        globs.tol = 1e-5;
        globs.maxIters = 200;

        Grid3D current(globs.nx, globs.ny, globs.nz);
        current.fill(0.0);
        current(globs.nx/2, globs.ny/2, globs.nz/2) = 100.0;

        LinearAlgebra lin(globs.maxIters);
        HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, lin);
        BoundaryConditions bc{{}, {}};

        for (int step = 0; step < 50; ++step) {
            Grid3D next = current;
            solver.step(current, next, globs, bc);
            current = next;
        }

        double center = current(globs.nx/2, globs.ny/2, globs.nz/2);
        double neighbor = current(globs.nx/2 + 1, globs.ny/2, globs.nz/2);

        if (center >= 100.0 || neighbor <= 0.0) {
            std::cout << "FAIL (center=" << center << ", neighbor=" << neighbor << ")\n";
            return 1;
        }
        std::cout << "PASS (center cooled to " << center << ")\n";
    }

    // Test 3: Rough analytical match (coarse grid)
    {
        std::cout << "Test 3: Analytical approximation... ";
        SimulationGlobals globs;
        globs.nx = globs.ny = globs.nz = 41;
        globs.dx = 0.025;
        globs.dt = 0.00005;
        globs.alpha = 0.01;
        globs.tol = 1e-4;
        globs.maxIters = 300;

        Grid3D current(globs.nx, globs.ny, globs.nz);
        current.fill(0.0);
        current(globs.nx/2, globs.ny/2, globs.nz/2) = 100.0;

        LinearAlgebra lin(globs.maxIters);
        HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, lin);
        BoundaryConditions bc{{}, {}};

        double t_sim = 0.0;
        for (int step = 0; step < 400; ++step) {
            Grid3D next = current;
            solver.step(current, next, globs, bc);
            current = next;
            t_sim += globs.dt;
        }

        double x = 5.0 * globs.dx;
        double expected = analytical_1d(x, t_sim, globs.alpha) * 100.0;
        double computed = current(globs.nx/2 + 5, globs.ny/2, globs.nz/2);

        if (std::abs(computed - expected) > 20.0) {
            std::cout << "FAIL (computed=" << computed << ", expected≈" << expected << ")\n";
            return 1;
        }
        std::cout << "PASS (difference within loose tolerance)\n";
    }
    //Test 4, 5, 6 Add, subtract, dot product two arrays and comapre CPU with GPU

    {   SimulationGlobals globs;
        std::size_t N = 20000;
        double min = -10000.0;
        double max = 10000.0;
        std::vector<double> randA =  randomVector(N, min, max);
        std::vector<double> randB =  randomVector(N, min, max);

        dim3 blockDim1D(globs.blockDimX*globs.blockDimY*globs.blockDimZ);
        dim3 gridDim1D((N+globs.blockDimX*globs.blockDimY*globs.blockDimZ-1)
                        /(globs.blockDimX*globs.blockDimY*globs.blockDimZ));

        
        double* devBVector = nullptr;
        double* devAVector = nullptr;
        double* devCVector = nullptr;
        double* devBlockSums = nullptr;

        CUDA_CHECK(cudaMalloc(&devAVector, N*sizeof(double)));
        CUDA_CHECK(cudaMalloc(&devBVector, N*sizeof(double)));
        CUDA_CHECK(cudaMalloc(&devCVector, N*sizeof(double)));
        CUDA_CHECK(cudaMalloc(&devBlockSums, gridDim1D.x*sizeof(double)));

        CUDA_CHECK(cudaMemcpy(devAVector, randA.data(), N*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(devBVector, randB.data(), N*sizeof(double), cudaMemcpyHostToDevice));

        std::cout << "Test 4: Add Two Arrays CPU and GPU match... ";
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
    }
    std::cout << "All tests pass! \n";
    return 0;
}
