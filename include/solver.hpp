#pragma once
#include "grid.hpp"
#include "simGlobals.hpp"
#include "boundaryConditions.hpp"
#include "linearAlgebra.hpp"
#include <iostream>

class HeatSolver {
    public:
        virtual ~HeatSolver() = default;
        
        virtual void step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs, const BoundaryConditions& bc) = 0;

        virtual const char* name() const = 0;

};