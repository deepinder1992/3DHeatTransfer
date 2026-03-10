#include "boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(std::array<BCType,6>types,std::array<double,6>values)
           :types_(types),values_(values){};


void BoundaryConditions::applyBCsToStencil(Grid3D& grid,double dx, double cond) const{
    const size_type nx = grid.nx();
    const size_type ny = grid.ny();
    const size_type nz = grid.nz();
    
    // -----------------------------
    // x faces
    // -----------------------------
    auto applyBC = [&](int face, double& cell, double interior, int sign ){ if (types_[face]==BCType::Dirichlet){cell=values_[face];}
            else if (types_[face]==BCType::Neumann){ cell= (sign*2*dx*values_[face]/cond)+interior;}};
    for (size_type k = 0; k < nz; ++k){
        for (size_type  j = 0; j < ny; ++j){
            //x min face
            applyBC(0, grid(0,j,k),grid(2,j,k),-1);
            //x max face
            applyBC(1, grid(nx-1,j,k),grid(nx-3,j,k),1);            
        }
    }
     // -----------------------------
    // y faces
    // -----------------------------
    for (size_type k = 0; k < nz; ++k)
        for (size_type i = 0; i < nx; ++i) {
            // ymin
            applyBC(2, grid(i,0,k), grid(i,2,k),-1);
            // ymax
            applyBC(3, grid(i,ny-1,k),grid(i,ny-3,k),1);
        }
    // -----------------------------
    // Z faces
    // -----------------------------
    for (size_type j = 0; j < ny; ++j)
        for (size_type i = 0; i < nx; ++i) {
            // zmin
            applyBC(4, grid(i,j,0),grid(i,j,2),-1);
            // zmax
            applyBC(5, grid(i,j,nz-1),grid(i,j,nz-3),1);
        }
}

void BoundaryConditions::applyBCsToRhsMatrix(size_type nx,
                                              size_type ny,
                                               size_type nz, 
                                                double dx,
                                                 double coeff,
                                                  double cond,
                                                   std::vector<double>& b) const
{   
 
    auto applyBC = [&](size_type faceInx, size_type row, int factor){
                    if(types_[faceInx]==BCType::Dirichlet)b[row] += 2.0*coeff*values_[faceInx];
                    else if(types_[faceInx]==BCType::Neumann)b[row] += (factor*2.0*coeff*dx/cond)*values_[faceInx]; };

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

void BoundaryConditions::applyBCsToStencilCUDA(double* grid, double dx, size_type nx, 
    size_type ny, size_type nz, double cond, dim3 gridCuda, dim3 blockCuda) const
    {         

        applyBCsToStencilKern<<<gridCuda, blockCuda>>>(grid, nx, ny, nz, dx, cond, types_, values_);
        cudaDeviceSynchronize();
    }