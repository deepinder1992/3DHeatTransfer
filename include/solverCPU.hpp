#pragma once
#include "solver.hpp"
#include "sparseMatrix.hpp"
#include "linearAlgebra.hpp"


class HeatSolverCPUStencil final : public HeatSolver {

    public:
        HeatSolverCPUStencil(double alpha, double dx, double dt, const LinearAlgebra& linAlgebra);

        void step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc) override;

        const char* name() const override {return "CPU Implicit Stencil";}

    private:
        double alpha_, dx_, dt_, coeff_;
        LinearAlgebra linAlgebra_;

};

class HeatSolverCPUMatrix final : public HeatSolver{

    public:
        HeatSolverCPUMatrix(const Grid3D& grid, std::size_t nx, std::size_t ny, std::size_t nz, double alpha, double dx, double dt, double k, const BoundaryConditions& bc, const LinearAlgebra& linAlgebra);
        
        void step(const Grid3D& current, Grid3D & next,const SimulationGlobals& globs, const BoundaryConditions& bc) override;
        const char* name() const override {return "CPU Impict Matrix";}
    
    private:
        SparseMatrix A_;
        double alpha_, dx_, dt_, coeff_, cond_;
        LinearAlgebra linAlgebra_;
};