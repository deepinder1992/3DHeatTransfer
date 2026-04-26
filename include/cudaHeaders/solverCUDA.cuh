#pragma once
#include "solver.hpp"
#include "sparseMatrix.hpp"
#include "cudaHeaders/boundaryConditions.cuh"
#include "cudaHeaders/linearAlgebra.cuh"
#include <iostream>

class HeatSolverCUDAStencil final: public HeatSolver{
    public:
        HeatSolverCUDAStencil(double alpha, double dx, double dt, const LinearAlgebraCUDA& linAlgebraCUDA);
        ~HeatSolverCUDAStencil()
                        {
                            if(devCurrent) cudaFree(devCurrent);
                            if(devNext) cudaFree(devNext);
                            if(devOld) cudaFree(devOld);
                            if(devMaxBlockError) cudaFree(devMaxBlockError);
                            if(devBcIndices) cudaFree(devBcIndices);
                            if(devIntIndices) cudaFree(devIntIndices);
                            if(devFaceTypes) cudaFree (devFaceTypes);
                            if(devNbrTypes) cudaFree(devNbrTypes); 
                            if(devNbrOffset) cudaFree(devNbrOffset);   
                        }

        void step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc) override;

        const char* name() const override {return "CUDA Implicit Stencil";}

    private:
        double alpha_, dx_,dt_,coeff_;
        LinearAlgebraCUDA linAlgebraCUDA_;
        BoundaryConditionsCUDA bcCUDA_;

        double *devCurrent = nullptr, *devNext = nullptr, *devOld = nullptr,
                *devMaxBlockError = nullptr;
        float  (*devCellNormals)[3] = nullptr;
        
        FaceType *devFaceTypes = nullptr;
        NeighbourType  *devNbrTypes = nullptr;

        std::size_t (*devIntIndices)[3] = nullptr, (*devBcIndices)[3] = nullptr, *devNbrOffset = nullptr;

        size_type devMemCurrGrdSize = 0, devMemNextGrdSize = 0, devMemOldGrdSize = 0,
                  devMemBcIndSize = 0, devMemIntIndSize = 0, devMemBlockErrorSize = 0,
                  devMemFaceTypeSize = 0, devMemNbrSize = 0,  devMemNbrOffsetSize =0,
                  devMemCellNormalSize = 0;
                 
};

class HeatSolverCUDAMatrix final: public HeatSolver{
    public:
        HeatSolverCUDAMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k,
                                            const BoundaryConditions& bc, const LinearAlgebraCUDA& linAlgebraCUDA);


        void step(const Grid3D& current, Grid3D & next,const SimulationGlobals& globs, const BoundaryConditions& bc) override;

        const char* name() const override {return "CUDA Implicit Matrix";}

    private:
        SparseMatrix A_;
        double alpha_, dx_, dt_, coeff_, cond_;
        LinearAlgebraCUDA linAlgebraCUDA_;

};

