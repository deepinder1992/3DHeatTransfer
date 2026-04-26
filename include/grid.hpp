#pragma once
#include<vector>
#include <map> 
#include<cassert>
#include<cstddef>
#include <cmath>
#include "simGlobals.hpp"


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
        Grid3D(std::size_t nx, std::size_t ny, std::size_t nz, double dx);
        
        double& operator()(std::size_t i, std::size_t j, std::size_t k);
        const double&  operator()(std::size_t i, std::size_t j, std::size_t k) const;

        CellType& cellType(std::size_t i, std::size_t j,std::size_t k);
        const CellType& cellType(std::size_t i, std::size_t j,std::size_t k) const;

        FaceType& faceType(std::size_t i, std::size_t j,std::size_t k);
        const FaceType& faceType(std::size_t i, std::size_t j,std::size_t k) const;
        const std::vector<FaceType>& faceTypeVect() const {return faceType_;}

        Vector& cellFaceNormal(std::size_t i, std::size_t j,std::size_t k);
        const std::vector<Vector>& cellFaceNormals() const;
        Vector cellFaceNormalized(std::size_t i, std::size_t j,std::size_t k) const;

        const std::vector<std::array<std::size_t,3>>& activeIndices() const;
        const std::vector<std::array<std::size_t,3>>& boundaryIndices()const;
        const std::vector<std::array<std::size_t,3>>& interiorIndices()const;
        const std::vector<NeighbourType>& flatNeigbourTypes() const;
        const std::vector<std::size_t>& offsetsNeighbourTypes()const;

        std::vector<NeighbourType> getSolidNeighbours(std::size_t i, std::size_t j, std::size_t k) const;
        std::vector<NeighbourType> findSolidNeighbours(std::size_t i, std::size_t j, std::size_t k) const;

        std::size_t nx() const noexcept {return nx_;}
        std::size_t ny() const noexcept {return ny_;}
        std::size_t nz() const noexcept {return nz_;}
        double dx() const noexcept {return dx_;}
        
        void adjustGrid(const double maxStlEdge);
        Vector gridCent();

        std::size_t size() const noexcept{return data_.size();}
        double* data() noexcept {return data_.data();} 
        const double* data() const noexcept{return data_.data();}

        std::size_t numBoundaryCells() const { return numBoundaryCells_; }

        std::size_t totalCellsInGeometry() const { return numActiveCells_;}

        void compactLookup();
        const std::vector<std::size_t>& compactLookup() const;
        
        void fill (double value);
        
        void detectBoundaries();
        void constructNeigbourMap(SolverType solver);

        void diagnostics() const;
        void assignNoneCells();

    private:
        std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept;
        
        std::size_t nx_, ny_,nz_;
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
        
        std::size_t numActiveCells_=0, numBoundaryCells_=0, numSolidCells_=0;

};