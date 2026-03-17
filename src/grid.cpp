#include "grid.hpp"


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

size_type Grid3D::index(size_type i, size_type j, size_type k) const noexcept
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

Vector& Grid3D::cellFaceNormal(size_type i, size_type j,size_type k){
     return boundaryNormal_[index(i,j,k)];}

const Vector& Grid3D::cellFaceNormal(size_type i, size_type j,size_type k) const{
     return boundaryNormal_[index(i,j,k)];}

FaceType& Grid3D::faceType(size_type i, size_type j,size_type k){
        return faceType_[index(i,j,k)];
    }

const FaceType& Grid3D::faceType(size_type i, size_type j,size_type k) const{
        return faceType_[index(i,j,k)];
    }

void Grid3D::detectBoundaries(){
    for (size_type k=0; k< nz_; ++k){
        for(size_type j=0; j<ny_; ++j){
            for(size_type i=0; i<nx_; ++i){
                if (cellType(i,j,k)==CellType::SOLID) continue;
                
                if (findSolidNeigbour(i,j,k).size() !=0){
                    cellType(i,j,k) = CellType::BOUNDARY;
                    boundaryIndices_.push_back({i,j,k});}
            }
        }
    }
}

const std::vector<std::array<std::size_t,3>>  Grid3D::findSolidNeigbour(std::size_t i, std::size_t j, std::size_t k) const{
    std::vector<std::array<std::size_t,3>> solidNeighbours;
    if (i > 0     && cellType(i-1,j,k) == CellType::SOLID) solidNeighbours.push_back({i-1,j,k});
    else if (i < nx_-1 && cellType(i+1,j,k) == CellType::SOLID)  solidNeighbours.push_back({i+1,j,k});
    else if (j > 0     && cellType(i,j-1,k) == CellType::SOLID)  solidNeighbours.push_back({i,j-1,k});
    else if (j < ny_-1 && cellType(i,j+1,k) == CellType::SOLID)  solidNeighbours.push_back({i,j+1,k});
    else if (k > 0     && cellType(i,j,k-1) == CellType::SOLID)  solidNeighbours.push_back({i,j,k-1});
    else if (k < nz_-1 && cellType(i,j,k+1) == CellType::SOLID)  solidNeighbours.push_back({i,j,k+1});
    return solidNeighbours;
}