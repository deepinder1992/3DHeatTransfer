#pragma once
#include "grid.hpp"
#include <array>
#include "simGlobals.hpp"
// #include "kernel.cuh"
// #include "cuda_runtime.h"

class BoundaryConditions{
    public:
        BoundaryConditions(std::array<BCType,3>types, std::array<double,3>values);
        
        std::array<BCType,3> types() const noexcept {return types_;}
        std::array<double,3> values() const noexcept {return values_;}
        
        void applyBCsToStencil(Grid3D& grid, double dx, double cond) const;
        
        void applyBCsToRhsMatrix(const Grid3D& grid, size_type nx, size_type ny, double dx,
                                     double cond, double coeff, std::vector<double>& b) const;

        // void applyBCsToStencilCUDA(double* grid, double* oldGrid, double dx, size_type nx, 
        //             size_type ny, size_type nz, size_type(*bcIndices)[3], FaceType* faceTypes, std::size_t nBcCells,
        //             NeighbourType* devNbrTypes, std::size_t* devNbrOffset, float (*devCellNormals)[3], double cond, dim3 gridCuda, dim3 blockCuda) const;
        
        int sign(Vector cellNormal, const std::array<size_type,3>&  cellIdxs, std::size_t i, std::size_t j, std::size_t k) const;
 
    private:
        std::array<BCType,3> types_;
        std::array<double,3> values_;
};