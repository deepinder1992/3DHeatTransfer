#include "grid.hpp"
#include <iostream>

Grid3D::Grid3D(size_type nx, size_type ny, size_type nz, double dx)
:nx_(nx),ny_(ny),nz_(nz), dx_(dx), data_(nx*ny*nz), cellType_(nx*ny*nz, CellType::INTERIOR),
faceType_(nx*ny*nz, FaceType::NONE),boundaryNormal_(nx*ny*nz)
    {
        assert(nx>0 && ny > 0 && nz >0);
    }

double& Grid3D::operator()(size_type i, size_type j, size_type k)
    {  assert(i<nx_ && j< ny_ && k <nz_);
        return data_[index(i,j,k)];
    }


const double& Grid3D::operator()(size_type i, size_type j, size_type k) const
    {  assert(i<nx_ && j< ny_ && k <nz_);
        return data_[index(i,j,k)];
    }


void Grid3D::fill(double value)
    {
        std::fill(data_.begin(),data_.end(),value);
    }

std::size_t Grid3D::index(size_type i, size_type j, size_type k) const noexcept
    {
        return i+ nx_*(j+ ny_*k);
    }

const CellType& Grid3D::cellType(size_type i, size_type j,size_type k) const{
        return cellType_[index(i,j,k)];
    }

CellType& Grid3D::cellType(size_type i, size_type j,size_type k){
        return cellType_[index(i,j,k)];
    }

const std::vector<std::array<std::size_t,3>>& Grid3D::boundaryIndices() const{
    return boundaryIndices_;
}

const std::vector<std::size_t>& Grid3D::interiorIdxs() const{ 
    return interiorIdxs_;
}

Vector& Grid3D::cellFaceNormal(size_type i, size_type j,size_type k){
     return boundaryNormal_[index(i,j,k)];}

Vector Grid3D::cellFaceNormalized(size_type i, size_type j,size_type k) const{
     return boundaryNormal_[index(i,j,k)].normalize();}

FaceType& Grid3D::faceType(size_type i, size_type j,size_type k){
        return faceType_[index(i,j,k)];
    }

const FaceType& Grid3D::faceType(size_type i, size_type j,size_type k) const{
        return faceType_[index(i,j,k)];
    }

void Grid3D::adjustGrid(const double maxStlEdge) {
    dx_ = maxStlEdge/(nx_);
 }

Vector Grid3D::gridCent(){
    //we have cube grid assumption in which the arbitrary shape lies
    //assume lx=ly=lz
    float cent = 0.5f*static_cast<float>(nx_)*static_cast<float>(dx_);
    return {cent, cent, cent};
 }

void Grid3D::detectBoundaries(){
    auto makeBoundary = [&](std::size_t i, std::size_t j, std::size_t k) {
        cellType(i,j,k) = CellType::BOUNDARY;
        boundaryIndices_.push_back({i,j,k});
        ++numBoundaryCells_;  };

    auto makeInterior = [&](std::size_t i, std::size_t j, std::size_t k) {
        interiorIdxs_.push_back(index(i,j,k));
        ++numInteriorCells_;  };

    for (size_type k=0; k< nz_; ++k){
        for(size_type j=0; j<ny_; ++j){
            for(size_type i=0; i<nx_; ++i){
                if (cellType(i,j,k)==CellType::SOLID) continue;
                    
                makeInterior(i,j,k);
                
                if (findSolidNeigbour(i,j,k).size() !=0){
                    makeBoundary(i,j,k);}                
            }
        }
    }

    std::cout <<"INNNNN" <<numInteriorCells_<< std::endl;
    std::cout <<"BNNNNN" <<numBoundaryCells_<< std::endl;
    std::cout <<"SNNNNN" <<numBoundaryCells_<< std::endl;
}

const std::vector<NeighbourType>  Grid3D::findSolidNeigbour(std::size_t i, std::size_t j, std::size_t k) const{
    std::vector<NeighbourType> solidNeighbours;

    if (i > 0     && cellType(i-1,j,k) == CellType::SOLID) solidNeighbours.push_back(NeighbourType::X_PREV);
    else if (i < nx_-1 && cellType(i+1,j,k) == CellType::SOLID)  solidNeighbours.push_back(NeighbourType::X_NEXT);
    else if (j > 0     && cellType(i,j-1,k) == CellType::SOLID)  solidNeighbours.push_back(NeighbourType::Y_PREV);
    else if (j < ny_-1 && cellType(i,j+1,k) == CellType::SOLID)  solidNeighbours.push_back(NeighbourType::Y_NEXT);
    else if (k > 0     && cellType(i,j,k-1) == CellType::SOLID)  solidNeighbours.push_back(NeighbourType::Z_PREV);
    else if (k < nz_-1 && cellType(i,j,k+1) == CellType::SOLID)  solidNeighbours.push_back(NeighbourType::Z_NEXT);
    else if (i==0)     solidNeighbours.push_back(NeighbourType::X_PREV); 
    else if (i==nx_-1) solidNeighbours.push_back(NeighbourType::X_NEXT);
    else if (j==0)     solidNeighbours.push_back(NeighbourType::Y_PREV);
    else if (j==ny_-1) solidNeighbours.push_back(NeighbourType::Y_NEXT);
    else if (k==0)     solidNeighbours.push_back(NeighbourType::Z_PREV);
    else if (k==nz_-1) solidNeighbours.push_back(NeighbourType::Z_NEXT);
    return solidNeighbours;
}

void Grid3D::diagnostics()const{
    size_type numInteriorCells=0, numBoundaryCells=0, numSolidCells=0,
              numInletCells = 0, numOutletCells=0, numWallCells=0,
              numNoneCells=0;
    for (size_type k=0; k< nz_; ++k){
        for(size_type j=0; j<ny_; ++j){
            for(size_type i=0; i<nx_; ++i){

                if(cellType(i,j,k) == CellType::BOUNDARY) ++numBoundaryCells;
                if(cellType(i,j,k) == CellType::SOLID) ++numSolidCells;
                if(cellType(i,j,k) == CellType::INTERIOR) ++numInteriorCells;
                if(faceType(i,j,k) == FaceType::INLET) ++numInletCells;
                if(faceType(i,j,k) == FaceType::OUTLET) ++numOutletCells;
                if(faceType(i,j,k) == FaceType::WALL) ++numWallCells;
                if(cellType(i,j,k) == CellType::BOUNDARY && faceType(i,j,k) == FaceType::NONE) ++numNoneCells;
                if(cellType(i,j,k) == CellType::SOLID){//&&(i==0||j==0||k==0||i==nx_-1||j==ny_-1||k==nz_-1)){
                    std::cout <<i<<" "<<j<<" "<<k<<std::endl;
                }
            }
        }
    }
    std::cout << "\n......Grid Stats......" << std::endl;
    std::cout << "Total Cells:      " << nx_ * ny_ * nz_ << std::endl;
    std::cout << "Interior Cells:   " << numInteriorCells << std::endl;
    std::cout << "Boundary Cells:   " << numBoundaryCells << std::endl;
    std::cout << "Solid Cells:      " << numSolidCells << std::endl;
    std::cout << "Inlet Cells:      " << numInletCells << std::endl;
    std::cout << "Outlet Cells:     " << numOutletCells << std::endl;
    std::cout << "Wall Cells:       " << numWallCells << std::endl;
    std::cout << "Unassigned Cells: " << numNoneCells <<"\n" <<std::endl;
}

void Grid3D::assignNoneCells(){
    for (size_type k=0; k< nz_; ++k){
        for(size_type j=0; j<ny_; ++j){
            for(size_type i=0; i<nx_; ++i){
                if(cellType(i,j,k) == CellType::BOUNDARY && faceType(i,j,k) == FaceType::NONE){
                 faceType(i,j,k)=FaceType::WALL;
                }
            }
        }
    }
}
