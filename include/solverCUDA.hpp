#pragma once
#include "solver.hpp"
#include "linearAlgebra.hpp"
#include "sparseMatrix.hpp"
#include <iostream>

class HeatSolverCUDAStencil final: public HeatSolver{
    public:
        HeatSolverCUDAStencil(double alpha, double dx, double dt, const LinearAlgebra& linAlgebra);
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


class HeatSolverCUDAMatrix final: public HeatSolver{
    public:
        HeatSolverCUDAMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k,
                                            const BoundaryConditions& bc, const LinearAlgebra& linAlgebra);


        void step(const Grid3D& current, Grid3D & next,const SimulationGlobals& globs, const BoundaryConditions& bc) override;

        const char* name() const override {return "CUDA Implicit Matrix";}

    private:
        SparseMatrix A_;
        double alpha_, dx_, dt_, coeff_, cond_;
        LinearAlgebra linAlgebra_;

};

