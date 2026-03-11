#include "kernel.cuh"
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <numeric> 
#include <vector>


__global__  void implicitJacobiKernel(double* oldVal, double* newVal, double* currentVal, std::size_t nx, std::size_t ny, std::size_t nz, double coeff_)
    {
        std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
        std::size_t k = blockIdx.z * blockDim.z + threadIdx.z;

        if(i>0 && i<nx-1 && j>0 && j<ny-1 && k>0 && k<nz-1 ){
            std::size_t idx = i + j*nx + k*nx*ny;
            double rhs = currentVal[idx];

            double sum = oldVal[(i-1) + j*nx +k*nx*ny] + oldVal[(i+1)+j*nx +k*nx*ny]
                       + oldVal[i + (j-1)*nx +k*nx*ny] + oldVal[i+(j+1)*nx +k*nx*ny]
                       + oldVal[i + j*nx +(k-1)*nx*ny] + oldVal[i+j*nx +(k+1)*nx*ny];
            
            newVal[idx] = (rhs + coeff_*sum)/(1+6*coeff_);        
        }

    }


__global__ void addSubtract(double* a , double* b , double* c, double fac, std::size_t N, double sign){
            std::size_t stride = blockDim.x*gridDim.x;
            std::size_t gTid = blockDim.x*blockIdx.x + threadIdx.x;
            for (std::size_t i = gTid; i< N; i+=stride) {c[i]= a[i]+sign*fac*b[i];}           
    }


__global__ void dotBlock(const double* a, const double* b, double* blockSum, std::size_t N)
{
    extern __shared__ double sData[];

    std::size_t tid = threadIdx.x;
    std::size_t gTid = threadIdx.x + blockIdx.x * blockDim.x;
    std::size_t stride = blockDim.x * gridDim.x;

    double temp = 0.0;

    // Each thread accumulates its stride
    for (std::size_t i = gTid; i < N; i += stride)
        temp += a[i] * b[i];

    sData[tid] = temp;
    __syncthreads();

    // Reduction in shared memory
    for (std::size_t s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
            sData[tid] += sData[tid + s];
        __syncthreads();
    }

    // First thread writes block sum
    if (tid == 0)
        blockSum[blockIdx.x] = sData[0];
}
 

__global__ void sparseMultiply (const double* values, const std::size_t* cols, const std::size_t* rowPtr, const double*x, double* y, std::size_t N){
            std::size_t gTid = threadIdx.x + blockDim.x*blockIdx.x;
            std::size_t stride = blockDim.x*gridDim.x;
            for(std::size_t i = gTid; i<N; i+=stride){
                double temp = 0.0;
                for (std::size_t idx = rowPtr[i]; idx <rowPtr[i+1];++idx){
                    temp+=values[idx]*x[cols[idx]];
                }
                y[i]= temp;
            }            
    }

__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError, std::size_t N, std::size_t nx, std::size_t ny){
        extern __shared__ double sData[];

        std::size_t tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z*blockDim.x*blockDim.y;
        std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
        std::size_t k = blockIdx.z * blockDim.z + threadIdx.z;
        std::size_t gTid = i + j*nx + k*nx*ny;


        if (gTid < N)sData[tid] = fabs(oldVal[gTid]-newVal[gTid]);
        else sData[tid] = 0.0;
        __syncthreads();

        for (std::size_t s = blockDim.x*blockDim.y*blockDim.z/2; s>0; s>>=1){
            if(tid<s) sData[tid] = fmax(sData[tid], sData[tid+s]);
            __syncthreads();
        }

        if (tid == 0) maxBlockError[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = sData[0];

}

__global__ void arraySumReduction (double* a, double* blockSum, std::size_t n){
            extern __shared__ double sData[];
            std::size_t gTid = threadIdx.x + blockDim.x*blockIdx.x;
            std::size_t tid = threadIdx.x;
            sData[tid] = (gTid<n)?a[gTid]:0;
            __syncthreads();
            for (std::size_t s = blockDim.x/2; s>0; s>>=1){
                if(tid<s) sData[tid]+=sData[tid+s];
                __syncthreads();
            }
            if(tid==0)blockSum[blockIdx.x] = sData[0];
        }  

double arraySum(double* a , std::size_t n){
        dim3 blockDim(256);   
        double* blockSum;
        dim3 gridDim((n+blockDim.x-1)/blockDim.x);
        cudaMalloc(&blockSum,gridDim.x*sizeof(double));
        std::size_t i = n;
        double* d_input = a;
        while(true){ 
            gridDim = dim3((i+blockDim.x-1)/blockDim.x);
            std::size_t sharedMemSize = blockDim.x*sizeof(double);
            arraySumReduction<<<gridDim,blockDim,sharedMemSize>>>(d_input,blockSum,i);
            d_input = blockSum;
            i = gridDim.x;
            if(i<128){break;}
        }
        std::vector<double> hostBlockSum(i);
        cudaMemcpy(hostBlockSum.data(), d_input, i*sizeof(double), cudaMemcpyDeviceToHost);
        cudaFree(blockSum);
        return std::accumulate(hostBlockSum.begin(),hostBlockSum.end(),0.0);

}


__global__ void applyBCsToStencilKern(double* grid, std::size_t nx, std::size_t ny, std::size_t nz, double dx, double cond,
                                     const BCType types_[6], const double values_[6]){
        
        std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
        std::size_t k = blockIdx.z * blockDim.z + threadIdx.z;
        
        if(i>=nx|| j>=ny|| k>=nz) return;
        
        auto applyBc = [=] __device__ (int face, double& cell, double interior, int sign){
            if (types_[face]==BCType::Dirichlet){cell = values_[face];}
                  else if (types_[face]==BCType::Neumann){ cell = (sign*2*dx*values_[face]/cond)+interior;} };
        //x min face
        if(i==0) applyBc(0, grid[0+j*nx+k*ny*nx], grid[2+j*nx+k*ny*nx], -1);
        //x max face
        else if(i==nx-1) applyBc(1, grid[nx-1+j*nx+k*ny*nx], grid[nx-3+j*nx+k*ny*nx], 1);
       
        //y min face
        else if(j==0) applyBc(2, grid[i+0*nx+k*ny*nx], grid[i+2*nx+k*ny*nx], -1);

        //y max face
        else if(j==ny-1) applyBc(3, grid[i+(ny-1)*nx+k*ny*nx], grid[i+(ny-3)*nx+k*ny*nx], 1);
        
        //z min face
        else if(k==0)applyBc(4, grid[i+j*nx+0*ny*nx], grid[i+j*nx+2*ny*nx], -1);
        //z max face
        else if(k==nz-1)applyBc(5, grid[i+j*nx+(nz-1)*ny*nx], grid[i+j*nx+(nz-3)*ny*nx], 1);
}
