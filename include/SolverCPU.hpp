#pragma once
#include "solver.hpp"

class HeatSolverCPU final : public HeatSolver {

    public:
        HeatSolverCPU(double alpha, double dx, double dt);

        void step(const Grid3D& current, Grid3D next) override;

        const char* name() const override {return "CPU";}

    private:
        double alpha_, dx_, dt_, coeff_;

}