#pragma once
#include "solver.hpp"
#include <cuda_runtime.h>

class HeatSolverCUDAStencil final: public HeatSolver{
    public:
        HeatSolverCUDAStencil(double alpha, double dx, double dt, const LinearAlgebra& linAlgebra):
                                alpha_(alpha), dx_(dx),dt_(dt),linAlgebra_(linAlgebra)
                            {  assert(alpha> 0.0);
                               assert(dx > 0.0);
                               assert (dt > 0.0);
                               coeff_ = alpha_*dt_/(dx_*dx_);
                        };
        ~HeatSolverCUDAStencil()
                        {
                            if(devCurrent) cudaFree(devCurrent);
                            if(devNext) cudaFree(devNext);
                            if(devOld) cudaFree(devOld);
                            if(devMaxBlockError) cudaFree(devMaxBlockError);
                        }

        void step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc) override;

        const char* name() const override {return "CUDA Implicit Stencil";}

    private:
        double alpha_, dx_,dt_,coeff_;
        LinearAlgebra linAlgebra_;
        double* devCurrent = nullptr;
        double* devNext = nullptr;
        double* devOld = nullptr;
        double* devMaxBlockError = nullptr;

        size_type devMemCurrGrdSize = 0;
        size_type devMemNextGrdSize = 0;
        size_type devMemOldGrdSize = 0;
        size_type devMemBlockErrorSize = 0;

};


// class HeatSolverCUDAMatrix final: public HeatSolver{
//     public:
//         HeatSolverCUDAMatrix(double alpha, double dx, double dt):alpha_(alpha), dx_(dx),dt_(dt)
//                             {  assert(alpha> 0.0);
//                                assert(dx > 0.0);
//                                assert (dt > 0.0);
//                                coeff_ = alpha_*dt_/(dx_*dx_);
//                         };
//         ~HeatSolverCUDAStencil();

//         void step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc) override;

//         const char* name() const override {return "CUDA Implicit Matrix";}

//     private:
//         double alpha_, dx_,dt_,coeff_;
//         double* devCurrent = nullptr;

//         size_type devMemCurrGrdSize = 0;

// };

