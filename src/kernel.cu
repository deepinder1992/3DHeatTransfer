#include "kernel.cuh"
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <numeric> 
#include <vector>


__global__  void implicitJacobiKernel(double* oldVal, double* newVal, double* currentVal, int nx, int ny, int nz, double coeff_)
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


__global__ void addSubtract(double* a , double* b , double* c, double fac, int N, double sign){
            std::size_t stride = blockDim.x*gridDim.x;
            std::size_t gTid = blockDim.x*blockIdx.x + threadIdx.x;
            for (std::size_t i = gTid; i< N; i+=stride) {c[i]= a[i]+sign*fac*b[i];}           
    }

// __global__ void dotBlock (const double* a, const double* b, double* blockSum, int N){
//             extern __shared__ double sData[];
//             std::size_t gTid = threadIdx.x + blockDim.x*blockIdx.x;
//             std::size_t tid = threadIdx.x;
//             double temp = (gTid<N)?a[gTid]*b[gTid]:0.0;
//             sData[tid] = temp;
//             for (std::size_t s = blockDim.x/2; s>0; s>>=1){
//                 if(tid<s) sData[tid]+=sData[tid+s];
//                 __syncthreads();
//             }
//             if(tid==0) blockSum[blockIdx.x] = sData[0];
//         } 
__global__ void dotBlock(const double* a, const double* b, double* blockSum, int N)
{
    extern __shared__ double sData[];

    int tid = threadIdx.x;
    int gTid = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;

    double temp = 0.0;

    // Each thread accumulates its stride
    for (int i = gTid; i < N; i += stride)
        temp += a[i] * b[i];

    sData[tid] = temp;
    __syncthreads();

    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
            sData[tid] += sData[tid + s];
        __syncthreads();
    }

    // First thread writes block sum
    if (tid == 0)
        blockSum[blockIdx.x] = sData[0];
}
 

// __global__ void arrayAtomicAdd (double* a, double* result, int n){
//             extern __shared__ double sData[];
//             std::size_t gTid = threadIdx.x + blockDim.x*blockIdx.x;
//             std::size_t tid = threadIdx.x;
//             std::size_t stride = gridDim.x*blockDim.x;
//             double temp = 0.0;
//             for (std::size_t i = gTid; i<n; i+=stride){
//                 temp+=a[i];
//             }
//             sData[tid] = temp;
//             __syncthreads();

//             for (std::size_t s = blockDim.x/2; s>0; s>>=1){
//                 if(tid<s) sData[tid]+=sData[tid+s];
//                 __syncthreads();
//             }
//             if(tid==0)atomicAdd(result,sData[0]);
// }   



__global__ void sparseMultiply (const double* values, const std::size_t* cols, const std::size_t* rowPtr, const double*x, double* y, int N){
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

        for (std::size_t s = blockDim.x*blockDim.y*blockDim.z/2; s>0; s>>=1){
            if(tid<s) sData[tid] = fmax(sData[tid], sData[tid+s]);
            __syncthreads();
        }

        if (tid == 0) maxBlockError[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = sData[0];

}

__global__ void arraySumReduction (double* a, double* blockSum, int n){
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

double arraySum(double* a , int n){
        dim3 blockDim(256);   
        double* blockSum;
        dim3 gridDim((n+blockDim.x-1)/blockDim.x);
        cudaMalloc(&blockSum,gridDim.x*sizeof(double));
        std::size_t i = n;
        double* d_input = a;
        while(true){ 
            gridDim = dim3((i+blockDim.x-1)/blockDim.x);
            long int sharedMemSize = blockDim.x*sizeof(double);
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
