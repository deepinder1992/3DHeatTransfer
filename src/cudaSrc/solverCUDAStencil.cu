#include "cudaHeaders/solverCUDA.cuh"
#include <stdexcept>
#include <cuda_runtime.h>
#include "cudaHeaders/kernel.cuh"

HeatSolverCUDAStencil::HeatSolverCUDAStencil(double alpha, double dx, double dt, const LinearAlgebraCUDA& linAlgebraCUDA):
                                alpha_(alpha), dx_(dx),dt_(dt),linAlgebraCUDA_(linAlgebraCUDA)
                            {  assert(alpha> 0.0);
                               assert(dx > 0.0);
                               assert (dt > 0.0);
                               coeff_ = alpha_*dt_/(dx_*dx_);
                        };

void HeatSolverCUDAStencil::step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc){

    std::size_t N =  current.size();
    const std::size_t nx = current.nx();
    const std::size_t ny = current.ny();
    const std::size_t nz = current.nz();
    const std::size_t nIntIdxs = (current.interiorIndices()).size();
    const std::size_t nBcIdxs = (current.boundaryIndices()).size();
    const std::size_t nSolidNbrIdxs = (current.flatNeigbourTypes()).size();

    ::allocateMemory(devCurrent, devMemCurrGrdSize, N);
    ::allocateMemory(devNext, devMemNextGrdSize, N);
    ::allocateMemory(devOld, devMemOldGrdSize, N);
    ::allocateMemory(devFaceTypes, devMemFaceTypeSize, N);
    ::allocateMemory(devCellNormals, devMemCellNormalSize, N);

    ::allocateMemory(devBcIndices, devMemBcIndSize, nBcIdxs);
    ::allocateMemory(devIntIndices, devMemIntIndSize, nIntIdxs);
    ::allocateMemory(devNbrTypes, devMemNbrSize, nSolidNbrIdxs );
    ::allocateMemory(devNbrOffset, devMemNbrOffsetSize, nBcIdxs+1);
   
    cudaMemcpy(devCurrent, current.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devOld, next.data(), N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devNext, next.data(), N*sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(devIntIndices, current.interiorIndices().data(), 
                nIntIdxs*sizeof(*devIntIndices), cudaMemcpyHostToDevice);

    cudaMemcpy(devBcIndices, current.boundaryIndices().data(), 
            nBcIdxs*sizeof(*devBcIndices), cudaMemcpyHostToDevice);

    cudaMemcpy(devFaceTypes, current.faceTypeVect().data(), 
            N*sizeof(FaceType), cudaMemcpyHostToDevice);

    cudaMemcpy(devCellNormals, current.cellFaceNormals().data(), 
            N*sizeof(*devCellNormals), cudaMemcpyHostToDevice);

    cudaMemcpy(devNbrTypes, current.flatNeigbourTypes().data(), 
        nSolidNbrIdxs*sizeof(NeighbourType), cudaMemcpyHostToDevice);

    cudaMemcpy(devNbrOffset, current.offsetsNeighbourTypes().data(), 
        (nBcIdxs+1)*sizeof(std::size_t), cudaMemcpyHostToDevice);
    
    dim3 blockDims(globs.blockDim);
    dim3 gridDims((nIntIdxs+globs.blockDim-1)/globs.blockDim);


    std::size_t numBlocks = gridDims.x;

    ::allocateMemory(devMaxBlockError, devMemBlockErrorSize, numBlocks);
    
    for (int iter = 0; iter<linAlgebraCUDA_.maxIters();++iter){
        BoundaryConditionsCUDA::applyBCsToStencilCUDA(bc.types().data(), bc.values().data(), devOld, devNext, dx_,nx, ny, nz, devBcIndices, devFaceTypes, nBcIdxs,
            devNbrTypes, devNbrOffset, devCellNormals, globs.k, gridDims, blockDims);   

        linAlgebraCUDA_.implicitJacobiCUDA(devOld, devNext, devCurrent, devIntIndices, nIntIdxs, nx, ny, nz, coeff_, gridDims, blockDims);

        if (globs.verbosity & SimulationGlobals::VERB_MEDIUM){
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "CUDA error Jacobi: " << cudaGetErrorString(err) << std::endl;
            }
        }
        
        std::size_t sharedMemSize = blockDims.x*sizeof(double);
        linAlgebraCUDA_.maxErrorCUDA(devOld, devNext, devMaxBlockError, devIntIndices, nIntIdxs, nx, ny, gridDims, blockDims, sharedMemSize);

        if (globs.verbosity & SimulationGlobals::VERB_MEDIUM){
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "CUDA error maxError Call: " << cudaGetErrorString(err) << std::endl;
            }
        }

        double hostMaxBlockError[numBlocks];
        cudaMemcpy(hostMaxBlockError, devMaxBlockError, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

        double maxErr = 0.0;
        for (std::size_t b = 0; b < numBlocks; ++b)
            maxErr = std::max(maxErr, hostMaxBlockError[b]);
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<maxErr<< std::endl;
            ++globs.totalIters;    
        }     
           
        if (maxErr<globs.tol)break;

        std::swap(devOld,devNext);
        linAlgebraCUDA_.adjustMaxItersIfNeeded(iter);
    }   
    cudaMemcpy(next.data(), devOld, N*sizeof(double), cudaMemcpyDeviceToHost);  
}


