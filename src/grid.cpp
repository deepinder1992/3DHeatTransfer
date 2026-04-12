#include "grid.hpp"
#include <algorithm>
#include <iostream>

Grid3D::Grid3D(size_type nx, size_type ny, size_type nz, double dx)
:nx_(nx),ny_(ny),nz_(nz), dx_(dx), data_(nx*ny*nz), cellType_(nx*ny*nz, CellType::INTERIOR),
faceType_(nx*ny*nz, FaceType::NONE),boundaryNormal_(nx*ny*nz), compactLookup_(nx*ny*nz, INVALID),
offsetsNbrTypes_(1, 0)
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

const std::vector<std::array<std::size_t,3>>& Grid3D::activeIndices() const{ 
    return activeIndices_;
    }

const std::vector<std::array<std::size_t,3>>& Grid3D::interiorIndices() const{ 
    return interiorIndices_;
    }

Vector& Grid3D::cellFaceNormal(size_type i, size_type j,size_type k){
     return boundaryNormal_[index(i,j,k)];
    }

const std::vector<Vector>& Grid3D::cellFaceNormals() const{
     return boundaryNormal_;
    }

Vector Grid3D::cellFaceNormalized(size_type i, size_type j,size_type k) const{
     return boundaryNormal_[index(i,j,k)].normalize();
    }

FaceType& Grid3D::faceType(size_type i, size_type j,size_type k){
        return faceType_[index(i,j,k)];
    }

const FaceType& Grid3D::faceType(size_type i, size_type j,size_type k) const{
        return faceType_[index(i,j,k)];
    }

const std::vector<NeighbourType>& Grid3D::flatNeigbourTypes() const{
        return flatNbrTypes_;
    }
const std::vector<std::size_t>& Grid3D::offsetsNeighbourTypes() const{
        return offsetsNbrTypes_;
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

void Grid3D::compactLookup() {
     std::size_t counter = 0;
     for (auto& cell:activeIndices_){
        auto[i,j,k] = cell;
        compactLookup_[index(i,j,k)] = counter;
        ++counter;
     }
 }

const std::vector<std::size_t>& Grid3D::compactLookup() const {
        return compactLookup_;
 }

void Grid3D::detectBoundaries(){
    auto makeBoundary = [&](std::size_t i, std::size_t j, std::size_t k) {
        cellType(i,j,k) = CellType::BOUNDARY;
        boundaryIndices_.push_back({i,j,k});
        ++numBoundaryCells_;  };

    auto makeActive = [&](std::size_t i, std::size_t j, std::size_t k) {
        activeIndices_.push_back({i,j,k});
        ++numActiveCells_;  };
    
    auto makeInterior = [&](std::size_t i, std::size_t j, std::size_t k) {
        interiorIndices_.push_back({i,j,k}); };

    for (size_type k=0; k< nz_; ++k){
        for(size_type j=0; j<ny_; ++j){
            for(size_type i=0; i<nx_; ++i){
                if (cellType(i,j,k)==CellType::SOLID) continue;
                    
                makeActive(i,j,k);
                
                if (findSolidNeighbours(i,j,k).size() !=0){
                    makeBoundary(i,j,k);}    
                else{
                    makeInterior(i,j,k);
                }            
            }
        }
    }
    //make look up internal cells 
    compactLookup();
}

void Grid3D::constructNeigbourMap(SolverType solver){
    std::size_t idx = INVALID;


    std::vector<NeighbourType> solidNebrVect;
    for (auto ijk:boundaryIndices_){
        auto [i,j,k] =ijk;

        idx = index(i,j,k);

        auto removeCondition = [&](const NeighbourType& neighbour){
            int  s  = static_cast<int>(neighbour);
            auto ic = i + interiorOffsets[s][0];
            auto jc = j + interiorOffsets[s][1];
            auto kc = k + interiorOffsets[s][2];
            return cellType(ic,jc,kc)==CellType::SOLID;};
          
        solidNebrVect =  findSolidNeighbours(i, j, k);
        solidNebrVect.erase(
            std::remove_if(solidNebrVect.begin(), solidNebrVect.end(),
                            [&](const NeighbourType n){ 
                                return removeCondition(n);}), 
            solidNebrVect.end());

        solidNebrMap_[idx] = solidNebrVect;
        if(solidNebrMap_[idx].size()==0){
            std::cout<<"Cell centered at "<<i<<", "<<j<<", "<<k
            <<" is boundary but no solid neigbours found!"<<std::endl;}

        if (solver==SolverType::CUDA_STENCIL){
            flatNbrTypes_.insert(flatNbrTypes_.end(),
                                    solidNebrMap_[idx].begin(),
                                    solidNebrMap_[idx].end());

            offsetsNbrTypes_.push_back(flatNbrTypes_.size());
        }
    }
    //just a sanity check
    assert(solidNebrMap_.size()==numBoundaryCells_);

}


std::vector<NeighbourType> Grid3D::findSolidNeighbours(std::size_t i, std::size_t j, std::size_t k) const{
        std::vector<NeighbourType> solidsVect;
        if (i > 0     && cellType(i-1,j,k) == CellType::SOLID)  solidsVect.push_back(NeighbourType::X_PREV);
        if (i < nx_-1 && cellType(i+1,j,k) == CellType::SOLID)  solidsVect.push_back(NeighbourType::X_NEXT);
        if (j > 0     && cellType(i,j-1,k) == CellType::SOLID)  solidsVect.push_back(NeighbourType::Y_PREV);
        if (j < ny_-1 && cellType(i,j+1,k) == CellType::SOLID)  solidsVect.push_back(NeighbourType::Y_NEXT);
        if (k > 0     && cellType(i,j,k-1) == CellType::SOLID)  solidsVect.push_back(NeighbourType::Z_PREV);
        if (k < nz_-1 && cellType(i,j,k+1) == CellType::SOLID)  solidsVect.push_back(NeighbourType::Z_NEXT);
        if (i==0)      solidsVect.push_back(NeighbourType::X_PREV); 
        if (i==nx_-1)  solidsVect.push_back(NeighbourType::X_NEXT);
        if (j==0)      solidsVect.push_back(NeighbourType::Y_PREV);
        if (j==ny_-1)  solidsVect.push_back(NeighbourType::Y_NEXT);
        if (k==0)      solidsVect.push_back(NeighbourType::Z_PREV);
        if (k==nz_-1)  solidsVect.push_back(NeighbourType::Z_NEXT);
        return solidsVect;
    }

std::vector<NeighbourType>  Grid3D::getSolidNeighbours(std::size_t i, std::size_t j, std::size_t k) const{
    auto it = solidNebrMap_.find(index(i,j,k));
    if (it == solidNebrMap_.end()) {
        throw std::runtime_error("No solid neigbour at index: "
            + std::to_string(i) +" "+ std::to_string(i)+ " " + std::to_string(k));
    }
    return it->second;
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


