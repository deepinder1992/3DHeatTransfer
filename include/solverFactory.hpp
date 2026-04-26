#include "solverCPU.hpp"

#ifdef ENABLE_CUDA
#include "cudaHeaders/solverCUDA.cuh"
#endif 
   
std::unique_ptr<HeatSolver> CreateSolver(SolverType solverType, double alpha, Grid3D& current, double dt,  double k, int maxIters, BoundaryConditions& bc){  
    switch (solverType) {
        case SolverType::CPU_STENCIL: {
            LinearAlgebra linAlgebra(maxIters);        
            return  std::make_unique<HeatSolverCPUStencil>(alpha, current.dx(), dt, linAlgebra);
        }
        case SolverType::CPU_MATRIX: {
            LinearAlgebra linAlgebra(maxIters);
            return std::make_unique<HeatSolverCPUMatrix>(current, current.nx(), current.ny(), current.nz(), alpha, current.dx(),
                                dt, k, bc, linAlgebra); }
        #ifdef ENABLE_CUDA
            case SolverType::CUDA_STENCIL: {
                LinearAlgebraCUDA linAlgebraCUDA(maxIters);
                retunr std::make_unique<HeatSolverCUDAStencil>(alpha, current.dx(), dt, linAlgebraCUDA);
            }
            case SolverType::CUDA_MATRIX: {
                LinearAlgebraCUDA linAlgebraCUDA(globs.maxIters);
                return std::make_unique<HeatSolverCUDAMatrix>(current, current.nx(), current.ny(), current.nz(), alpha, current.dx(),
                                dt, k, bc, linAlgebra); 
                }
        #endif 
        default:
            throw std::runtime_error("Unknown solver selected! If you selected solver 3, 4 make sure to build with CUDA");
        }
    }