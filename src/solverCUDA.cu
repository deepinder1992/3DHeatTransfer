#include "solverCUDA.hpp"
#include <stdexcept>


template<typename T>
inline void allocateMemory(T*& ptr, size_type& allocatedSize, size_type N){
    if (N != allocatedSize){
        if (ptr) cudaFree(ptr);
        if (cudaMalloc(&ptr, N*sizeof(T))!= cudaSuccess)
            throw std::runtime_error("cudaMalloc Failed");
        allocatedSize = N;
    }
}

__global__  void implicitJacobiKernel(double* oldVal, double* newVal, double* currentVal, int nx, int ny, int nz, double coeff_)
    {
        std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
        std::size_t k = blockIdx.z * blockDim.z + threadIdx.z;

        if(i>0 && i<nx-1 && j>0 && j<ny-1 && k>0 && k<nz-1 ){
            size_t idx = i + j*nx + k*nx*ny;
            double rhs = currentVal[idx];

            double sum = oldVal[(i-1) + j*nx +k*nx*ny] + oldVal[(i+1)+j*nx +k*nx*ny]
                       + oldVal[i + (j-1)*nx +k*nx*ny] + oldVal[i+(j+1)*nx +k*nx*ny]
                       + oldVal[i + j*nx +(k-1)*nx*ny] + oldVal[i+j*nx +(k+1)*nx*ny];
            
            newVal[idx] = (rhs + coeff_*sum)/(1+6*coeff_);        
        }

    }

__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError, int N, int nx, int ny){
        extern __shared__ double sData[];

        std::size_t tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z*blockDim.x*blockDim.y;
        std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
        std::size_t k = blockIdx.z * blockDim.z + threadIdx.z;
        std::size_t gTid = i + j*nx + k*nx*ny;


        if (gTid < N)sData[tid] = fabs(oldVal[gTid]-newVal[gTid]);
        else sData[tid] = 0.0;
        __syncthreads();

        for (unsigned int s = blockDim.x*blockDim.y*blockDim.z/2; s>0; s>>=1){
            if(tid<s) sData[tid] = fmax(sData[tid], sData[tid+s]);
            __syncthreads();
        }

        if (tid == 0) maxBlockError[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = sData[0];

}


void HeatSolverCUDAStencil::step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc){

    size_type N = current.size();
    const std::size_t nx = current.nx();
    const std::size_t ny = current.ny();
    const std::size_t nz = current.nz();

    ::allocateMemory(devCurrent, devMemCurrGrdSize, N);
    ::allocateMemory(devNext, devMemNextGrdSize, N);
    ::allocateMemory(devOld, devMemOldGrdSize, N);
   
    cudaMemcpy(devCurrent, current.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devOld, next.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devNext, next.data(), N*sizeof(double), cudaMemcpyHostToDevice);

    
    dim3 blockDims(globs.blockDimX,globs.blockDimY, globs.blockDimZ);
    dim3 gridDims((nx+globs.blockDimX-1)/globs.blockDimX,
                     (ny+globs.blockDimY-1)/globs.blockDimY,
                        (nz+globs.blockDimZ-1)/globs.blockDimZ);

    size_type numBlocks = gridDims.x*gridDims.y*gridDims.z;

    ::allocateMemory(devMaxBlockError, devMemBlockErrorSize, numBlocks);
    
    for (int iter = 0; iter<=globs.maxIters;++iter){
        implicitJacobiKernel<<<gridDims, blockDims>>>(devOld, devNext, devCurrent, nx, ny, nz, coeff_);
        cudaDeviceSynchronize();
        
        if (globs.verbosity & SimulationGlobals::VERB_MEDIUM){
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "CUDA error Jacobi: " << cudaGetErrorString(err) << std::endl;
            }
        }
        long int sharedMemSize = blockDims.x*blockDims.y*blockDims.z*sizeof(double);
        maxError<<<gridDims, blockDims,sharedMemSize>>>(devOld, devNext, devMaxBlockError, N, nx, ny);
        if (globs.verbosity & SimulationGlobals::VERB_MEDIUM){
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "CUDA error maxError Call: " << cudaGetErrorString(err) << std::endl;
            }
        }

        double hostMaxBlockError[numBlocks];
        cudaMemcpy(hostMaxBlockError, devMaxBlockError, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

        double maxErr = 0.0;
        for (int b = 0; b < numBlocks; ++b)
            maxErr = std::max(maxErr, hostMaxBlockError[b]);
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<maxErr<< std::endl;
            ++globs.totalIters;    
        }     
           
        if (maxErr<globs.tol)break;

        std::swap(devOld,devNext);
    }   
    cudaMemcpy(next.data(), devOld, N*sizeof(double), cudaMemcpyDeviceToHost);

    bc.applyBCsToStencil(next, globs.dx, globs.k);   

    
}


