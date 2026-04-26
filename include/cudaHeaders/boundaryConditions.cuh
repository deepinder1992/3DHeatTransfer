#pragma once
#include "boundaryConditions.hpp"
#include "kernel.cuh"
#include "cuda_runtime.h"

class BoundaryConditionsCUDA{
    public:
        static void applyBCsToStencilCUDA(const BCType* types, const double* values, double* grid, double* oldGrid, double dx, std::size_t nx, 
                    std::size_t ny, std::size_t nz, std::size_t(*bcIndices)[3], FaceType* faceTypes, std::size_t nBcCells,
                    NeighbourType* devNbrTypes, std::size_t* devNbrOffset, float (*devCellNormals)[3], double cond, dim3 gridCuda, dim3 blockCuda);
};