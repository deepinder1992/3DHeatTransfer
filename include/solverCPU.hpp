#pragma once
#include "solver.hpp"
#include "sparseMatrix.hpp"


class HeatSolverCPUStencil final : public HeatSolver {

    public:
        HeatSolverCPUStencil(double alpha, double dx, double dt);

        void step(const Grid3D& current, Grid3D& next, const SimulationGlobals& globs, const BoundaryConditions& bc) override;

        const char* name() const override {return "CPU Implicit Stencil";}

    private:
        double alpha_, dx_, dt_, coeff_;

};

class HeatSolverCPUMatrix final : public HeatSolver{

    public:
        HeatSolverCPUMatrix(size_type nx, size_type ny, size_type nz, double alpha, double dx, double dt, double k, const BoundaryConditions& bc);
        
        void step(const Grid3D& current, Grid3D & next,const SimulationGlobals& globs, const BoundaryConditions& bc) override;
        const char* name() const override {return "CPU Impict Matrix";}
    
    private:
        SparseMatrix A_;
        double alpha_, dx_, dt_, coeff_, cond_;

};