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


__global__ void subtract(double*a , double* b , double* c, int N);


__global__ void add(double*a , double* b , double* c, int N);


__global__ void dotBlock (double* a, double* b, double* blockSum, int N);


__global__ void sparseMultiply (const double* values, const std::size_t* cols, const std::size_t* rowPtr, const double*x, double* y, int N);


__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError, int N, int nx, int ny);







// void conjugateGradientCUDA(const SparseMatrix& A,
//                        const std::vector<double>& b,
//                        std::vector<double>& x,
//                          const SimulationGlobals& globs)
// {
//     int N = b.size();
//     std::vector<double> r(N), p(N), Ap(N);

//     sparseMultiply<<<globs.gridDim1D, globs.blockDim1D>>>(A, A.colIndex(), A.rowPtr(), x.data(), Ap.data(), N);
//     for (int i = 0; i < N; ++i)
//         r[i] = b[i] - Ap[i];

//     p = r;
//     double rsold = dot(r, r);
//         for (int iter = 0; iter < globs.maxIters; ++iter)
//     {   
//         sparseMultiply<<<globs.gridDim1D, globs.blockDim1D>>>(A, A.colIndex(), A.rowPtr(), p.data(), Ap.data(), N);

//         double alpha = rsold / dot(p, Ap);

//         for (int i = 0; i < N; ++i){
//             x[i] += alpha * p[i];
//             r[i] -= alpha * Ap[i];}

//         double rsnew = dot(r, r);
//         double sqrtRsnew = std::sqrt(rsnew);
//         if (globs.verbosity & SimulationGlobals::VERB_HIGH){
//             std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<sqrtRsnew << std::endl;
//             ++globs.totalIters;}
//         if (sqrtRsnew < globs.tol) break;

//         for (int i = 0; i < N; ++i)
//             p[i] = r[i] + (rsnew / rsold) * p[i];

//         rsold = rsnew;
//     }
// }
// void launchdotBlock(double* a, double* b, double* blockSum, int N)
// {
//     dotBlock (a, b, blockSum, N);
// }
