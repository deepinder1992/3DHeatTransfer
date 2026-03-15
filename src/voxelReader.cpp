#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "voxelReader.hpp"

bool rayIntersectsTriangle(const Vector& origin, const Vector& direction, const Triangle& tri) {
    constexpr double epsilon = 1e-6;
    constexpr double epsilon2 = 1e-5;
    Vector normal = tri.normal;
    double nDenom = normal.dot(direction);
    if(abs(nDenom) < epsilon){return false;}// with in parallel threshold
    double t = normal.dot(tri.v0-origin)/(normal.dot(direction));
    if(t <= epsilon2){return false;}// with

    Vector p = origin + direction*t;

    Vector e1 = tri.v1 - tri.v0;
    Vector e2 = tri.v2 - tri.v0;
    Vector vp = p - tri.v0;

    double dot00 = e2.dot(e2);
    double dot01 = e2.dot(e1);
    double dot02 = e2.dot(vp);
    double dot11 = e1.dot(e1);
    double dot12 = e1.dot(vp);

    double bCentDenom = (dot00 * dot11 - dot01 * dot01);
    if (abs(bCentDenom) <= epsilon){return false;} //degenrate triangle

    double u = (dot11 * dot02 - dot01 * dot12) / bCentDenom;
    double v = (dot00 * dot12 - dot01 * dot02) / bCentDenom;

    // inside triangle test
    if (u >= -epsilon && v >= -epsilon && (u + v) <= 1.0+epsilon)
        return true;

    return false;
}

std::string addSuffix(const std::string& filename, 
                                        const std::string& suffix) {
    size_t dotPosition = filename.find_last_of('.');
    std::string base = filename.substr(0, dotPosition);
    std::string ext  = filename.substr(dotPosition);    
    return base +'_'+ suffix + ext;
}

double distanceFromCentroid(double x, double y, double z, const Triangle& tri)
{
    double centX = (tri.v0.x + tri.v1.x + tri.v2.x)/3.0;
    double centY = (tri.v0.y + tri.v1.y + tri.v2.y)/3.0;
    double centZ = (tri.v0.z + tri.v1.z + tri.v2.z)/3.0;
    double distx = centX-x;
    double disty = centY-y;
    double distz = centZ-z;
    return std::sqrt(distx*distx + disty*disty + distz*distz);
}


VoxelReader::VoxelReader(const std::string& fileName, Grid3D& grid){
        std::vector<Triangle> triangles;
        static_assert(sizeof(Vector) == 3 * sizeof(float), "Vec3 must have no padding");
        std::cout << "Reading Stl file..\n";
        loadBinaryStl(fileName, triangles);

        std::cout <<"Started Voxelizing..\n";

        std::cout << "Ray tracing..\n";
        voxelizeGrid(grid, triangles);
        
        std::cout << "Detecting boundaries..\n";
        grid.detectBoundaries();

        //add inlet,outlet and wall to file name
        std::string inletFile = addSuffix(fileName,"inlet");
        std::string outletFile = addSuffix(fileName,"outlet");
        std::string wallFile = addSuffix(fileName,"wall");
        
        std::cout << "Applying inlet patch..\n";
        std::vector<Triangle> inletTriangles;
        loadBinaryStl(inletFile, inletTriangles);
        voxelizePatch(grid, inletTriangles, FaceType::INLET);

        std::cout << "Applying outlet patch..\n";
        std::vector<Triangle> outletTriangles;
        loadBinaryStl(outletFile, outletTriangles);
        voxelizePatch(grid, outletTriangles, FaceType::OUTLET);
        
        std::cout << "Applying wall patch..\n";
        std::vector<Triangle> wallTriangles;
        loadBinaryStl(wallFile, wallTriangles);
        voxelizePatch(grid, wallTriangles, FaceType::WALL);
        std::cout <<"Finished Voxelizing!\n";
}

bool VoxelReader::loadBinaryStl(const std::string& fileName, std::vector<Triangle>& triangles ){
    if(!voxelReadBinaryStl(fileName, triangles)){
        std::cerr << "Error: failed to read binary STL file " << fileName << "\n";
        return false;}
    return true;}

bool VoxelReader::voxelReadBinaryStl(const std::string& fileName, std::vector<Triangle>& triangles){
    std::ifstream fin(fileName,std::ios::binary);
    if(!fin.is_open()) return false;
    char header[80];
    fin.read(header,80);
    uint32_t numTriangles;
    triangles.resize(numTriangles);
    fin.read(reinterpret_cast<char*>(&numTriangles), sizeof(uint32_t));
    for (uint32_t t=0; t<numTriangles;++t){
        Triangle& tri = triangles[t];
        fin.read(reinterpret_cast<char*>(&tri.normal), 3*sizeof(float));
        fin.read(reinterpret_cast<char*>(&tri.v0), 3*sizeof(float));
        fin.read(reinterpret_cast<char*>(&tri.v1), 3*sizeof(float));
        fin.read(reinterpret_cast<char*>(&tri.v2), 3*sizeof(float));
        char attr[2];
        fin.read(attr,2);
    }
    fin.close();
    return true;
}

void VoxelReader::voxelizeGrid(Grid3D& grid, const std::vector<Triangle>& triangles){
    std::size_t nx = grid.nx();
    std::size_t ny = grid.ny();
    std::size_t nz = grid.nz();
    //uniform grid assumption todo make it non uniform
    double dx = grid.dx();
    double dy = grid.dx();
    double dz = grid.dx();

    auto isInside = [&](const Vector& voxelCenter, const std::vector<Triangle>& triangles){
        Vector rayDir = {1.0,0,0};
        int crossCount = 0;
        for (const Triangle& tri:triangles){
            if(rayIntersectsTriangle(voxelCenter,rayDir, tri)){
                crossCount++;}
         } 
        return (crossCount%2)==1; 
    };

    for (std::size_t k = 0; k<nz; ++k){
        for (std::size_t j=0; j<ny; ++j){    
            for (std::size_t i=0; i<nx; ++i){
                float x = (i+0.5)*dx;
                float y = (j+0.5)*dy;
                float z = (k+0.5)*dz;
                Vector origin = {x,y,z};
                grid.cellType(i,j,k) = (isInside(origin, triangles))?CellType::INTERIOR:CellType::SOLID;                    
            }
        }
    }
}

void VoxelReader::voxelizePatch(Grid3D& grid, const std::vector<Triangle>& triangles,const FaceType faceType){
    std::size_t nx = grid.nx();
    std::size_t ny = grid.ny();
    std::size_t nz = grid.nz();
    //uniform grid assumption todo make it non uniform
    double dx = grid.dx();
    double dy = grid.dx();
    double dz = grid.dx();

    for(const Triangle& tri:triangles){
        float xmin = std::min(tri.v0.x,std::min(tri.v1.x,tri.v2.x));
        float xmax = std::max(tri.v0.x,std::max(tri.v1.x,tri.v2.x));
        float ymin = std::min(tri.v0.y,std::min(tri.v1.y,tri.v2.y));
        float ymax = std::max(tri.v0.y,std::max(tri.v1.y,tri.v2.y));
        float zmin = std::min(tri.v0.z,std::min(tri.v1.z,tri.v2.z));
        float zmax = std::max(tri.v0.z,std::max(tri.v1.z,tri.v2.z));
        std::size_t minI =  static_cast<std::size_t>(floor(xmin/dx));
        std::size_t maxI = std::min(static_cast<std::size_t>(floor(xmax/dx)),nx-1);
        std::size_t minJ =  static_cast<std::size_t>(floor(ymin/dy));
        std::size_t maxJ = std::min(static_cast<std::size_t>(floor(ymax/dy)),ny-1);
        std::size_t minK =  static_cast<std::size_t>(floor(zmin/dz));
        std::size_t maxK = std::min(static_cast<std::size_t>(floor(zmax/dz)),nz-1);
        for (std::size_t k = minK; k<=maxK; ++k){
            for (std::size_t j=minJ; j<maxJ; ++j){    
                for (std::size_t i=minI; i<maxI; ++i){
                    double x = (i+0.5)*dz;
                    double y = (j+0.5)*dy;
                    double z = (k+0.5)*dz;
                    if (grid.cellType(i,j,k) == CellType::BOUNDARY && distanceFromCentroid(x, y, z, tri)<dx){
                        {
                            grid.faceType(i,j,k) = faceType;
                        } 
                    }                 
                }
            }
        }}
}


