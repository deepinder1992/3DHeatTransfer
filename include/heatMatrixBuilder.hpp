#pragma once 
#include "sparseMatrix.hpp"
#include <array>

inline SparseMatrix implicitMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double coeff, const BoundaryConditions& bc){
        size_type N = grid.totalCellsInGeometry();
        SparseMatrix A(N);

        auto& values = A.values();
        auto& cols   = A.colIndex();
        auto& rowPtr = A.rowPtr();

        const std::vector<std::size_t>& compactLookup = grid.compactLookup();

        rowPtr[0] =0;
        auto addNeighbor = [&](size_type ni,size_type nj, size_type nk,double val){
                size_type nidx = ni+nx*(nj+ny*nk);
                size_type col = compactLookup[nidx];
                if (col!= INVALID){
                    values.push_back(val);
                    cols.push_back(col);
                }
        };       

        auto applyBC = [&](size_type faceInx,size_type ni, size_type nj, size_type nk, double& diagVal){
                if (bc.types()[faceInx]==BCType::Dirichlet)diagVal+=coeff;
                else if (bc.types()[faceInx]==BCType::Neumann)addNeighbor(ni,nj,nk,-coeff);
                else diagVal+=coeff;
        };         

        for (auto& cell:grid.activeIndices()){
            auto [i,j,k] = cell;
            size_type idx =  i + nx*(j+ny*k);
            size_type row =  compactLookup[idx];

            double diag = 1+6.0*coeff;
            //handle boundaries
            if(grid.cellType(i,j,k) == CellType::BOUNDARY){
                const std::vector<NeighbourType> solidNeighbors = grid.getSolidNeighbours(i, j, k);
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