#include "cudaHeaders/kernel.cuh"
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <numeric> 
#include <vector>

__device__ __constant__ int interiorOffsetsGPU[6][3] =  { {+2,0,0} ,   {-3,0,0}, 
                                                      {0,+2,0} ,   {0,-3,0}, 
                                                      {0,0,+2} ,   {0,0,-3}};


__global__  void implicitJacobiKernel(double* oldVal, double* newVal, double* currentVal, std::size_t (*intIndices)[3], std::size_t nIntIdxs,
                                    std::size_t nx, std::size_t ny, std::size_t nz, double coeff_)
    {
        std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= nIntIdxs) return ;

        const auto& ijk = intIndices[tid];
        std::size_t i = ijk[0];
        std::size_t j = ijk[1];
        std::size_t k = ijk[2];

        std::size_t idx = i + j*nx + k*nx*ny;
        double rhs = currentVal[idx];

        double sum = oldVal[(i-1) + j*nx +k*nx*ny] + oldVal[(i+1)+j*nx +k*nx*ny]
                    + oldVal[i + (j-1)*nx +k*nx*ny] + oldVal[i+(j+1)*nx +k*nx*ny]
                    + oldVal[i + j*nx +(k-1)*nx*ny] + oldVal[i+j*nx +(k+1)*nx*ny];
        
        newVal[idx] = (rhs + coeff_*sum)/(1+6*coeff_);

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

__global__ void maxError( double* oldVal, double* newVal, double* maxBlockError,
                         std::size_t (*intIndices)[3], std::size_t nIntIdxs,
                         std::size_t nx, std::size_t ny){

    extern __shared__ double sData[];
    std::size_t tid  = blockIdx.x * blockDim.x + threadIdx.x;
    std::size_t lTid = threadIdx.x;

    if (tid >= nIntIdxs) sData[lTid] = 0.0;
    else{ const auto& ijk = intIndices[tid];
          std::size_t i = ijk[0];
          std::size_t j = ijk[1];
          std::size_t k = ijk[2];
          std::size_t idx = i + j*nx + k*nx*ny;
          sData[lTid] = fabs(oldVal[idx] - newVal[idx]);}

    __syncthreads();

    for (std::size_t s = blockDim.x / 2; s > 0; s >>= 1){
        if (lTid < s) sData[lTid] = fmax(sData[lTid], sData[lTid + s]);
        __syncthreads();
    }

    if (lTid == 0) maxBlockError[blockIdx.x] = sData[0];
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


__global__ void applyBCsToStencilKern(double* grid, double* oldGrid, std::size_t nx, std::size_t ny, std::size_t nz, double dx, 
                                    std::size_t (*bcIndices)[3], FaceType* faceTypes,std::size_t nBcCells,
                                    NeighbourType* nbrType, std::size_t* nbrOffset, float (*devCellNormals)[3],
                                    double  cond, const BCType types_[3], const double values_[3]){
        
        auto sign = [=] __device__ (float cellNormal[3],
                                    std::size_t i, std::size_t j, std::size_t k,
                                    std::size_t ic, std::size_t jc, std::size_t kc) {
                float radial[3] =  {static_cast<float>(i) - static_cast<float>(ic),
                                        static_cast<float>(j) - static_cast<float>(jc),
                                        static_cast<float>(k) - static_cast<float>(kc)};
                //dot product 
                if((radial[0]*cellNormal[0]+radial[1]*cellNormal[1]+radial[2]*cellNormal[2])>0) return 1;
                return -1; };
                
        std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

        if (tid>=nBcCells) return;
        const auto& ijk = bcIndices[tid];
        std::size_t i = ijk[0];
        std::size_t j = ijk[1];
        std::size_t k = ijk[2];
        std::size_t inx = i+j*nx+k*ny*nx;
        int faceNum = static_cast<int>(faceTypes[inx])-1;

        std::size_t start = nbrOffset[tid];
        std::size_t end   = nbrOffset[tid + 1];
        double weightBc = (start!=end)? 1.0/(static_cast<double>(end)-static_cast<double>(start)):1.0;
        
        //__syncthreads();
        grid[inx] = 0.0;            
        for (std::size_t idx = start; idx < end; ++idx) {
            if(types_[faceNum]==BCType::Dirichlet){
                grid[inx] +=weightBc*values_[faceNum];}
            else if (types_[faceNum]==BCType::Neumann){
                NeighbourType& neighbour = nbrType[idx];
                int  s  = static_cast<int>(neighbour);
                auto ic = i + interiorOffsetsGPU[s][0];
                auto jc = j + interiorOffsetsGPU[s][1];
                auto kc = k + interiorOffsetsGPU[s][2];

                int sign_ = sign(devCellNormals[inx], i,j,k,ic,jc,kc );
                grid[inx] += weightBc*((sign_*2*dx*values_[faceNum]/cond)+oldGrid[ic+jc*nx+kc*ny*nx]);}    
        }
}
