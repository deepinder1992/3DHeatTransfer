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

// __global__ void applyBcsStencil(Grid3D grid, std::size_t nx, std::size_t ny, std::size_t nz, double dx, double cond,
//                                     std::array<BCType,6> types_, std::array<double,6> values_){
        
//         std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
//         std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
//         std::size_t k = blockIdx.z * blockDim.z + threadIdx.z;
        
//         if(i>=nx|| j>=ny|| k>=nz) return;
        
//         //x min face
//         if(i==0){ if (types_[0]==BCType::Dirichlet){grid(0,j,k)=values_[0];}
//                   else if (types_[0]==BCType::Neumann){ grid(0,j,k)= (-2*dx*values_[0]/cond)+grid(2,j,k);}}
//         //x max face
//         if(i==nx-1){ if (types_[1]==BCType::Dirichlet){grid(nx-1,j,k)=values_[1];}
//                      else if (types_[1]==BCType::Neumann){ grid(nx-1,j,k)= (2*dx*values_[1]/cond)+grid(nx-3,j,k);} }
       
//         //y min face
//         if(j==0){ if (types_[2] == BCType::Dirichlet) grid(i,0,k) = values_[2];
//                   else if (types_[2] == BCType::Neumann) grid(i,0,k) =  (-2*dx*values_[2]/cond)+grid(i,2,k);}

//         //y max face
//         if(j==ny-1){ if (types_[3] == BCType::Dirichlet) grid(i,ny-1,k) = values_[3];
//                      else if (types_[3] == BCType::Neumann) grid(i,ny-1,k) =  (2*dx*values_[3]/cond)+grid(i,ny-3,k); }
        
//         //z min face
//         if(k==0){ if (types_[4] == BCType::Dirichlet) grid(i,j,0) = values_[4];
//             else if (types_[4] == BCType::Neumann) grid(i,j,0) =  (-2*dx*values_[4]/cond)+grid(i,j,2);}
//         //z max face
//         if(k==nz-1){ if (types_[5] == BCType::Dirichlet) grid(i,j,nz-1) = values_[5];
//                      else if (types_[5] == BCType::Neumann) grid(i,j,nz-1) =  (2*dx*values_[5]/cond)+grid(i,j,nz-3);}
//         return;
// }