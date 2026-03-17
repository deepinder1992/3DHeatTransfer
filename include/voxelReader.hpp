#pragma once
#include "grid.hpp"
#include <string>
#include <vector>


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
        static bool voxelReadBinaryStl(const std::string& fileName, std::vector<Triangle>& triangles);
        static void voxelizeGrid(Grid3D& grid, const std::vector<Triangle>& triangles);
        static void voxelizePatch(Grid3D& grid, const std::vector<Triangle>& triangles,const FaceType faceType);
        static void boundingBox(const Triangle& tri, float& minX,float& minY, float& minZ,
                                float& maxX,float& maxY, float& maxZ );
        static void shiftTriangles(std::vector<Triangle>& triangles);
    private:
        static float bBoxMinX, bBoxMaxX, bBoxMinY,bBoxMaxY, bBoxMinZ, bBoxMaxZ;    
};

