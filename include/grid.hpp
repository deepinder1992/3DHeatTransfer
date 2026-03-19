#pragma once
#include<vector>
#include<cassert>
#include<cstddef>
/////////////////////
////////////////
//////////
////////
///////
/////////
///////
#include <iostream>
////////////////

////////////////
///////////
///////////////////////
////////////////
//////////

///////
#include "cuda_runtime.h"
#include "simGlobals.hpp"


using size_type = std::size_t;

#pragma pack(push, 1)
struct Vector{
    float x, y, z;
    Vector operator- (const Vector& vectB)const {return {x-vectB.x, y-vectB.y, z-vectB.z};}
    Vector operator+ (const Vector& vectB)const {return {x+vectB.x, y+vectB.y, z+vectB.z};}
    Vector operator* (float s)const {return {x*s,y*s,z*s};}

    float dot(const Vector& vectB){return {x*vectB.x+ y*vectB.y+ z*vectB.z};}
};
#pragma pack(pop)


class Grid3D{
    public:
        Grid3D(size_type nx, size_type ny, size_type nz, double dx);
        
        double& operator()(size_type i, size_type j, size_type k);
        const double&  operator()(size_type i, size_type j, size_type k) const;

        CellType& cellType(size_type i, size_type j,size_type k);
        const CellType& cellType(size_type i, size_type j,size_type k) const;

        FaceType& faceType(size_type i, size_type j,size_type k);
        const FaceType& faceType(size_type i, size_type j,size_type k) const;

        Vector& cellFaceNormal(size_type i, size_type j,size_type k);
        const Vector& cellFaceNormal(size_type i, size_type j,size_type k) const;
    
        const std::vector<std::array<std::size_t,3>>& boundaryIndices()const;
        const  std::vector<std::array<std::size_t,3>> findSolidNeigbour(std::size_t i, std::size_t j, std::size_t k) const;

        size_type nx() const noexcept {return nx_;}
        size_type ny() const noexcept {return ny_;}
        size_type nz() const noexcept {return nz_;}
        double dx() const noexcept {return dx_;}
        
        void adjustGrid(const double maxStlEdge);
        Vector gridCent();

        size_type size() const noexcept{return data_.size();}
        double* data() noexcept {return data_.data();} 
        const double* data() const noexcept{return data_.data();}

        size_type numInteriorCells() const {return numInteriorCells_;}

        size_type numBoundaryCells() const { return numBoundaryCells_; }
        
        void fill (double value);
        
        void detectBoundaries();
    private:
        size_type index(size_type i, size_type j, size_type k) const noexcept;
        
        size_type nx_, ny_,nz_;
        double dx_;

        std::vector<double> data_;

        std::vector<CellType> cellType_;

        std::vector<FaceType> faceType_;
        
        std::vector<std::array<std::size_t,3>> boundaryIndices_;

        std::vector<Vector> boundaryNormal_;

        size_type numInteriorCells_=0, numBoundaryCells_=0;

};