#pragma once 
#include "sparseMatrix.hpp"
#include <array>

inline SparseMatrix implicitMatrix(const Grid3D& grid, std::size_t nx, std::size_t ny, std::size_t nz, double coeff, const BoundaryConditions& bc){
        std::size_t N = grid.totalCellsInGeometry();
        SparseMatrix A(N);

        auto& values = A.values();
        auto& cols   = A.colIndex();
        auto& rowPtr = A.rowPtr();

        const std::vector<std::size_t>& compactLookup = grid.compactLookup();

        rowPtr[0] =0;
        auto addNeighbor = [&](std::size_t ni,std::size_t nj, std::size_t nk,double val){
                std::size_t nidx = ni+nx*(nj+ny*nk);
                std::size_t col = compactLookup[nidx];
                if (col!= INVALID){
                    values.push_back(val);
                    cols.push_back(col);
                }
        };       

        auto applyBC = [&](std::size_t faceInx,std::size_t ni, std::size_t nj, std::size_t nk, double& diagVal){
                if (bc.types()[faceInx]==BCType::Dirichlet)diagVal+=coeff;
                else if (bc.types()[faceInx]==BCType::Neumann)addNeighbor(ni,nj,nk,-coeff);
                else diagVal+=coeff;
        };         

        for (auto& cell:grid.activeIndices()){
            auto [i,j,k] = cell;
            std::size_t idx =  i + nx*(j+ny*k);
            std::size_t row =  compactLookup[idx];

            double diag = 1+6.0*coeff;
            //handle boundaries
            if(grid.cellType(i,j,k) == CellType::BOUNDARY){
                const std::vector<NeighbourType> solidNeighbors = grid.findSolidNeighbours(i, j, k);
                int faceNum = static_cast<int>(grid.faceType(i,j,k))-1;
                for (NeighbourType neighbour :solidNeighbors){
                        if(neighbour==NeighbourType::X_PREV) applyBC(faceNum, i+1,j,k, diag);
                        if(neighbour==NeighbourType::X_NEXT) applyBC(faceNum,i-1,j,k,diag); 
                        if(neighbour==NeighbourType::Y_PREV) applyBC(faceNum,i,j+1,k,diag);        
                        if(neighbour==NeighbourType::Y_NEXT) applyBC(faceNum,i,j-1,k,diag);        
                        if(neighbour==NeighbourType::Z_PREV) applyBC(faceNum,i,j,k+1,diag);       
                        if(neighbour==NeighbourType::Z_NEXT) applyBC(faceNum,i,j,k-1,diag);}
                }


            values.push_back(diag);
            cols.push_back(row);
            
            if(i> 0)    addNeighbor(i-1,j,k,-coeff);
            if(i<nx-1)  addNeighbor(i+1,j,k,-coeff);
            if(j>0)     addNeighbor(i,j-1,k,-coeff);
            if(j<ny-1)  addNeighbor(i,j+1,k,-coeff);
            if(k>0)     addNeighbor(i,j,k-1,-coeff);
            if(k<nz-1)  addNeighbor(i,j,k+1,-coeff);

            rowPtr[row+1] = values.size();
        }

    return A;
}