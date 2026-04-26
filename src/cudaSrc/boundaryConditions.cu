#include "cudaHeaders/boundaryConditions.cuh"

void BoundaryConditionsCUDA::applyBCsToStencilCUDA(const BCType* types, const double* values, double* grid, double* oldGrid, double dx, std::size_t nx, 
    std::size_t ny, std::size_t nz, std::size_t(*bcIndices)[3], FaceType* faceTypes, std::size_t nBcCells,
    NeighbourType* devNbrTypes, std::size_t* devNbrOffset, float (*devCellNormals)[3], double cond, dim3 gridCuda, dim3 blockCuda)  
    {    
        // int numFaces = 3;
        // BCType types[numFaces];
        // for (int i = 0; i < numFaces; ++i)
        //     types[i] = types_[i];

        // double values[numFaces];
        // for (int i = 0; i < numFaces; ++i)
        //     values[i] = values_[i];

        applyBCsToStencilKern<<<gridCuda, blockCuda>>>(grid, oldGrid, nx, ny, nz, dx, bcIndices, faceTypes,
             nBcCells, devNbrTypes, devNbrOffset, devCellNormals, cond, types, values);
        cudaDeviceSynchronize();
    }