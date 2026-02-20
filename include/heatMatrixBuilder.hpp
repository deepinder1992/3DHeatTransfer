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
        auto applyBC = [&](size_type faceInx,size_type ni, size_type nj, size_type nk, double& diagVal){
                if (globs.types[faceInx]==BCType::Dirichlet)diagVal+=coeff;
                else if (globs.types[faceInx]==BCType::Neumann)addNeighbor(ni,nj,nk,-coeff);
                else diagVal+=coeff;
        };         

        for (size_type k =0 ; k< nz; ++k)
        for (size_type j =0 ; j< ny; ++j)
        for (size_type i =0 ; i< nx; ++i){
            size_type row =  i + nx*(j+ny*k);
            
            double diag = 1+6.0*coeff;
            
            if(i==0) applyBC(0, i+1,j,k, diag);// ----- X- -----
            if(i==nx-1) applyBC(1,i-1,j,k,diag); // ----- X+ -----
            if(j==0) applyBC(2,i,j+1,k,diag); // ----- X- -----       
            if(j==ny-1) applyBC(3,i,j-1,k,diag);// ----- Y+ -----        
            if(k==0)applyBC(4,i,j,k+1,diag);// ----- Z- -----       
            if(k==nz-1)applyBC(5,i,j,k-1,diag);// ----- Z+ -----


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