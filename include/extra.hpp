#pragma once 
#include "sparseMatrix.hpp"
#include <array>
#include "bcType.hpp"

void applyBoundaryConditions(size_type nx,
                             size_type ny,
                             size_type nz,
                             const SimulationGlobals& globs,
                             std::vector<double>& b)
{   
    double coeff = globs.alpha*globs.dt/(globs.dx*globs.dx);
    for (size_type k=0;k<nz;++k)
    for (size_type j=0;j<ny;++j)
    for (size_type i=0;i<nx;++i)
    {
        size_type row = i + nx*(j+ny*k);

        // ----- X- -----
        if(i==0)
        {
            if(globs.types[0]==BCType::Dirichlet)
                b[row] += 2.0*coeff*globs.values[0];

            else if(globs.types[0]==BCType::Neumann)
                b[row] -= (2.0*coeff*globs.dx/globs.k)*globs.values[0];
        }

        // ----- X+ -----
        if(i==nx-1)
        {
            if(globs.types[1]==BCType::Dirichlet)
                b[row] += 2.0*coeff*globs.values[1];

            else if(globs.types[1]==BCType::Neumann)
                b[row] += (2.0*coeff*globs.dx/globs.k)*globs.values[1];
        }

        // ----- Y- -----
        if(j==0)
        {
            if(globs.types[2]==BCType::Dirichlet)
                b[row] += 2.0*coeff*globs.values[2];

            else if(globs.types[2]==BCType::Neumann)
                b[row] -= (2.0*coeff*globs.dx/globs.k)*globs.values[2];
        }

        // ----- Y+ -----
        if(j==ny-1)
        {
            if(globs.types[3]==BCType::Dirichlet)
                b[row] += 2.0*coeff*globs.values[3];

            else if(globs.types[3]==BCType::Neumann)
                b[row] += (2.0*coeff*globs.dx/globs.k)*globs.values[3];
        }

        // ----- Z- -----
        if(k==0)
        {
            if(globs.types[4]==BCType::Dirichlet)
                b[row] += 2.0*coeff*globs.values[4];

            else if(globs.types[4]==BCType::Neumann)
                b[row] -= (2.0*coeff*globs.dx/globs.k)*globs.values[4];
        }

        // ----- Z+ -----
        if(k==nz-1)
        {
            if(globs.types[5]==BCType::Dirichlet)
                b[row] += 2.0*coeff*globs.values[5];

            else if(globs.types[5]==BCType::Neumann)
                b[row] += (2.0*coeff*globs.dx/globs.k)*globs.values[5];
        }
    }
}
