#include "boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(std::array<BCType,6>types,std::array<double,6>values)
           :types_(types),values_(values){};


void BoundaryConditions::apply(Grid3D& grid, int t) const{
    const size_type nx = grid.nx();
    const size_type ny = grid.ny();
    const size_type nz = grid.nz();


    //x faces 0 and 1
    for (size_type k = 0; k < nz; ++k){
        for (size_type  j = 0; j < ny; ++j){
            //x min face
            //apply dirchlet only once at t=0
            if (t== 0 && types_[0]==BCType::Dirichlet){grid(0,j,k)=values_[0];}
            else if (types_[0]==BCType::Neumann){ grid(0,j,k)=grid(1,j,k);}

            //x max face
            if (t== 0 && types_[1]==BCType::Dirichlet){grid(nx-1,j,k)=values_[1];}
            else if (types_[1]==BCType::Neumann){ grid(nx-1,j,k)=grid(nx-2,j,k);}
            
        }
    }
    for (size_type k = 0; k < nz; ++k)
        for (size_type i = 0; i < nx; ++i) {
            // ymin
            if (t== 0 && types_[2] == BCType::Dirichlet) grid(i,0,k) = values_[2];
            else if (types_[2] == BCType::Neumann) grid(i,0,k) = grid(i,1,k);

            // ymax
            if (t== 0 && types_[3] == BCType::Dirichlet) grid(i,ny-1,k) = values_[3];
            else if (types_[3] == BCType::Neumann) grid(i,ny-1,k) = grid(i,ny-2,k);
        }

    // -----------------------------
    // Z faces
    // -----------------------------
    for (size_type j = 0; j < ny; ++j)
        for (size_type i = 0; i < nx; ++i) {
            // zmin
            if (t== 0 && types_[4] == BCType::Dirichlet) grid(i,j,0) = values_[4];
            else if (types_[4] == BCType::Neumann) grid(i,j,0) = grid(i,j,1);

            // zmax
            if (t== 0 && types_[5] == BCType::Dirichlet) grid(i,j,nz-1) = values_[5];
            else if (types_[5] == BCType::Neumann) grid(i,j,nz-1) = grid(i,j,nz-2);
        }



}