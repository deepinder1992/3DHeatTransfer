#pragma once
#include "grid.hpp"
#include <string>
#include <vector>



#pragma pack(push, 1)
struct Vector{
    float x, y, z;
    Vector operator- (const Vector& vectB)const {return {x-vectB.x, y-vectB.y, z-vectB.z};}
    Vector operator+ (const Vector& vectB)const {return {x+vectB.x, y+vectB.y, z+vectB.z};}
    Vector operator* (float s)const {return {x*s,y*s,z*s};}

    float dot(const Vector& vectB){return {x*vectB.x+ y*vectB.y+ z*vectB.z};}
};
#pragma pack(pop)

struct Triangle{
    Vector normal;
    Vector v0;
    Vector v1;
    Vector v2;
};

class VoxelReader{
    public:
        VoxelReader(const std::string& fileName, Grid3D& grid);
        static bool loadBinaryStl(const std::string& fileName,std::vector<Triangle>& triangles);
    private:
        static bool voxelReadBinaryStl(const std::string& fileName, std::vector<Triangle>& triangles);
        static void voxelizeGrid(Grid3D& grid, const std::vector<Triangle>& triangles);
        static void voxelizePatch(Grid3D& grid, const std::vector<Triangle>& triangles,const FaceType faceType);
    
};

