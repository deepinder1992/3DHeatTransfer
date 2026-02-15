#pragma once
#include "grid.hpp"
#include "simGlobals.hpp"

class HeatSolver {
    public:
        virtual ~HeatSolver() = default;
        
        virtual void step(const Grid3D& current, Grid3D& next,const SimulationGlobals& globs) = 0;

        virtual const char* name() const = 0;

};