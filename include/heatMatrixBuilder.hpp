#pragma once 
#include "sparseMatrix.hpp"


SparseMatrix implicitMatrix( size_type nx, size_type ny, size_type nz, double coeff){
        size_type N = nx*ny*nz;
        SparseMatrix A(N);

        auto& values = A.values();
        auto& cols   = A.colIndex();
        auto& rowPtr = A.rowPtr();

        rowPtr[0] =0;
        
        for (int k =0 ; k< nz; ++k)
        for (int j =0 ; j< ny; ++j)
        for (int i =0 ; i< nx; ++i){
            int  row =  i + nx*(j+ny*k);
            
            values.push_back(1+6.0*coeff);
            cols.push_back(row);

            auto addNeighbor = [&](int ni,int nj, int nk){
                int col = ni+nx*(nj+ny*nk);
                values.push_back(-coeff);
                cols.push_back(col);
            };

            if(i> 0)    addNeighbor(i-1,j,k);
            if(i<nx-1)  addNeighbor(i+1,j,k);
            if(j>0)     addNeighbor(i,j-1,k);
            if(j<ny-1)  addNeighbor(i,j+1,k);
            if(k>0)     addNeighbor(i,j,k-1);
            if(k<nz-1)  addNeighbor(i,j,k+1);

            rowPtr[row+1] = values.size();

        }
    return A;
}