#pragma once
#include<vector>
#include <map> 
#include<cassert>
#include<cstddef>
#include <cmath>
#include "cuda_runtime.h"
#include "simGlobals.hpp"


using size_type = std::size_t;

#pragma pack(push, 1)
struct Vector{
    float x, y, z;

    Vector operator- (const Vector& vectB)const {return {x-vectB.x, y-vectB.y, z-vectB.z};}

    Vector operator+ (const Vector& vectB)const {return {x+vectB.x, y+vectB.y, z+vectB.z};}

    Vector operator* (float s)const {return {x*s,y*s,z*s};}

    Vector operator^ (const Vector& vectB)const {return { y * vectB.z - z * vectB.y,
                                                                 z * vectB.x - x * vectB.z,
                                                                 x * vectB.y - y * vectB.x };}


    float dot(const Vector& vectB){return {x*vectB.x+ y*vectB.y+ z*vectB.z};}
 
    float mag() const {return sqrtf(x*x + y*y + z*z);}

    Vector normalize() const {float m = mag();
        if (m<1e-8f)return{0.0,0.0,0.0};
        return {x/m,y/m,z/m};}
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
        const std::vector<FaceType>& faceTypeVect() const {return faceType_;}

        Vector& cellFaceNormal(size_type i, size_type j,size_type k);
        const std::vector<Vector>& cellFaceNormals() const;
        Vector cellFaceNormalized(size_type i, size_type j,size_type k) const;

        const std::vector<std::array<std::size_t,3>>& activeIndices() const;
        const std::vector<std::array<std::size_t,3>>& boundaryIndices()const;
        const std::vector<std::array<std::size_t,3>>& interiorIndices()const;
        const std::vector<NeighbourType>& flatNeigbourTypes() const;
        const std::vector<std::size_t>& offsetsNeighbourTypes()const;

        std::vector<NeighbourType> getSolidNeighbours(std::size_t i, std::size_t j, std::size_t k) const;
        std::vector<NeighbourType> findSolidNeighbours(std::size_t i, std::size_t j, std::size_t k) const;

        size_type nx() const noexcept {return nx_;}
        size_type ny() const noexcept {return ny_;}
        size_type nz() const noexcept {return nz_;}
        double dx() const noexcept {return dx_;}
        
        void adjustGrid(const double maxStlEdge);
        Vector gridCent();

        size_type size() const noexcept{return data_.size();}
        double* data() noexcept {return data_.data();} 
        const double* data() const noexcept{return data_.data();}

        size_type numBoundaryCells() const { return numBoundaryCells_; }

        size_type totalCellsInGeometry() const { return numActiveCells_;}

        void compactLookup();
        const std::vector<std::size_t>& compactLookup() const;
        
        void fill (double value);
        
        void detectBoundaries();
        void constructNeigbourMap(SolverType solver);

        void diagnostics() const;
        void assignNoneCells();

    private:
        std::size_t index(size_type i, size_type j, size_type k) const noexcept;
        
        size_type nx_, ny_,nz_;
        double dx_;

        std::vector<double> data_;

        std::vector<CellType> cellType_;
        std::vector<FaceType> faceType_;

        std::vector<std::array<std::size_t,3>> boundaryIndices_, activeIndices_, interiorIndices_;
        std::vector<Vector> boundaryNormal_;

        std::vector<std::size_t> compactLookup_;

        std::map<std::size_t, std::vector<NeighbourType>> solidNebrMap_;
        std::vector<NeighbourType> flatNbrTypes_;
        std::vector<std::size_t> offsetsNbrTypes_;
        
        size_type numActiveCells_=0, numBoundaryCells_=0, numSolidCells_=0;

};