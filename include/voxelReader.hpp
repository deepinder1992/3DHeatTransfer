#pragma once
#include "grid.hpp"
#include <string>
#include <vector>


struct Triangle{
    float normal[3];
    float v0[3];
    float v1[3];
    float v2[3];
}
struct Vector{
    double x, y,z;
    Vector operator- (const& vectB){x-vectB.x, y-vectB.y, z-vectC.z}; 
    Vector operator+ (const& vectB){x+vectB.x, y+vectB.y, z+vectC.z};
    Vector operator* (double s){x*s,y*s,z*s};

    double dot(const& vectB){x*vectB.x+y*vectB.y+z*vectB.z};
}
class VoxelReader{
    public:
        VoxelReader(const std::string& fileName, Grid3D& grid);
        static bool loadBinaryStl(const std::string& fileName, Grid3D& grid,std::vector<Triangle>& triangles);
        static bool rayIntersectsTriangle(const Vector& origin, const Vector& direction, const Triangle& tri);
        static double distanceFromCentroid(double x, double y, double z, const Triangle& tri);

    private:
        static bool readBinaryStl(const std::string fileName, std::vector<Triangle>& triangles);
        static void voxelizeGrid(Grid3D& grid, const std::vector<Triangle>& triangles)
    
};

