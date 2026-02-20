#pragma once 
#include "sparseMatrix.hpp"
#include <array>
#include "bcType.hpp"

void applyBoundaryConditionsRHSMatrix(size_type nx,
                             size_type ny,
                             size_type nz,
                             const SimulationGlobals& globs,
                             std::vector<double>& b)
{   
    double coeff = globs.alpha*globs.dt/(globs.dx*globs.dx);
    auto applyBC = [&](size_type faceInx, size_type row, int factor){
                    if(globs.types[faceInx]==BCType::Dirichlet)b[row] += 2.0*coeff*globs.values[faceInx];
                    else if(globs.types[faceInx]==BCType::Neumann)b[row] += (factor*2.0*coeff*globs.dx/globs.k)*globs.values[faceInx]; };

    // for (size_type k :{size_type(0), nz-1})
    // for (size_type j :{size_type(0), ny-1})
    // for (size_type i :{size_type(0), nx-1})
    for (size_type k =0 ; k< nz; ++k)
    for (size_type j =0 ; j< ny; ++j)
    for (size_type i =0 ; i< nx; ++i)
    {
        size_type row = i + nx*(j+ny*k);

        if(i==0) applyBC(0, row, -1);// ----- X- -----
        if(i==nx-1) applyBC(1,row,1); // ----- X+ -----
        if(j==0) applyBC(2,row,-1); // ----- X- -----       
        if(j==ny-1) applyBC(3,row,1);// ----- Y+ -----        
        if(k==0)applyBC(4,row,-1);// ----- Z- -----       
        if(k==nz-1)applyBC(5,row,1);// ----- Z+ -----

    }
}
