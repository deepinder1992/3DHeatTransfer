#pragma once
#include "grid.hpp"


class HeatSolver {
    public:
        virtual ~HeatSolver() = default;
        
        virtual void step(const Grid3D& current, Grid3D& next) = 0;

        virtual const char* name() const = 0;

}