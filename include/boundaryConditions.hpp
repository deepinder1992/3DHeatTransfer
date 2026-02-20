#pragma once
#include "grid.hpp"
#include "bcType.hpp"
#include <array>

class BoundaryConditions{
    public:
        BoundaryConditions(std::array<BCType,6>types, std::array<double,6>values);
        
        std::array<BCType,6> types() const noexcept {return types_;}
        std::array<double,6> values() const noexcept {return values_;}
        
        void applyBCsToStencil(Grid3D& grid, double dx, double cond) const;
        void applyBCsToRhsMatrix(size_type nx, size_type ny, size_type nz, double dx, double cond, double coeff, std::vector<double>& b) const;
        
    private:
        std::array<BCType,6> types_;
        std::array<double,6> values_;

};