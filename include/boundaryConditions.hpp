#pragma once
#include "grid.hpp"
#include "bcType.hpp"
#include <array>

class BoundaryConditions{
    public:
        BoundaryConditions(std::array<BCType,6>types,std::array<double,6>values);
          //  :types_(types),values_(values);
        void apply(Grid3D& grid) const;
    private:
        std::array<BCType,6> types_;
        std::array<double,6> values_;

};