#pragma once 
#include "sparseMatrix.hpp"
#include <array>
#include "bcType.hpp"

SparseMatrix implicitMatrix( size_type nx, size_type ny, size_type nz,const SimulationGlobals& globs){
        size_type N = nx*ny*nz;
        SparseMatrix A(N);

        auto& values = A.values();
        auto& cols   = A.colIndex();
        auto& rowPtr = A.rowPtr();

        double coeff = globs.alpha*globs.dt/(globs.dx*globs.dx);

        rowPtr[0] =0;
        auto addNeighbor = [&](size_type ni,size_type nj, size_type nk,double val){
                size_type col = ni+nx*(nj+ny*nk);
                values.push_back(val);
                cols.push_back(col);
            };


        for (size_type k =0 ; k< nz; ++k)
        for (size_type j =0 ; j< ny; ++j)
        for (size_type i =0 ; i< nx; ++i){
            size_type row =  i + nx*(j+ny*k);
            double diag = 1+6.0*coeff;
            if(i==0){
                if (globs.types[0]==BCType::Dirichlet)diag+=coeff;
                else if (globs.types[0]==BCType::Neumann)addNeighbor(i+1,j,k,-coeff);
                else diag+=coeff;

            }
            if(i==nx-1){
                if (globs.types[1]==BCType::Dirichlet)diag+=coeff;
                else if (globs.types[1]==BCType::Neumann)addNeighbor(i-1,j,k,-coeff);
                else diag+=coeff;
            }
            if(j==0){
                if (globs.types[2]==BCType::Dirichlet)diag+=coeff;
                else if (globs.types[2]==BCType::Neumann)addNeighbor(i,j+1,k,-coeff);
                else diag+=coeff;
            }
            if(j==ny-1){
                if (globs.types[3]==BCType::Dirichlet)diag+=coeff;
                else if (globs.types[3]==BCType::Neumann)addNeighbor(i,j-1,k,-coeff);
                else diag+=coeff;
            }
            if(k==0){
                if (globs.types[4]==BCType::Dirichlet)diag+=coeff;
                else if (globs.types[4]==BCType::Neumann)addNeighbor(i,j,k+1,-coeff);
                else diag+=coeff;
            }
            if(k==nz-1){
                if (globs.types[5]==BCType::Dirichlet)diag+=coeff;
                else if (globs.types[5]==BCType::Neumann)addNeighbor(i,j,k-1,-coeff);
                else diag+=coeff;
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