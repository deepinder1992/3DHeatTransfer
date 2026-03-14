#include <fstream>
#include <iostream>
#include <algorithm>
#include "voxelReader.hpp"

VoxelReader::VoxelReader(const std::string& fileName, Grid3D& grid){
        loadBinaryStl(fileName, grid, triagles )
        voxelizeGrid(Grid3D& grid, const triangles);
        grid.detectBoundaries();
        //add inlet,outlet and wall to file name
        std::string inletFile = addSuffix(fileName,'inlet')
        std::string outletFile = addSuffix(fileName,'outlet')
        std::string wallFile = addSuffix(fileName,'wall')
        
        loadBinaryStl(inletFile, grid, triangles)

        loadBinaryStl(outletFile, grid, triangles)

        loadBinaryStl(wallFile, grid, triangles)
}

static bool VoxelReader::loadBinaryStl(const std::string& fileName, Grid3D& grid,std::vector<Triangle>& triangles ) {
    std::vector<Triangle> triangles;
    if(!readBinaryStl(const std::string fileName, triangles))
    {std::cerr << "Error: failed to read binary STL file " << filename << "\n";
    return false;}
    return triangles;

}

static bool VoxelReader::voxelReadBinaryStl(const std::string fileName, std::vector<Triangle>& triagles){
    std::ifstream fin(filename,std::ios::binary)
    if(!fin.is_open()) return false;
    char header[80];
    fin.read(header,80);
    uint32_t numTriangles;
    fin.read(reinterpret_cast<char*>&numTriangles, sizeof(uint32_t));
    for (uint32_t t; t<numTriangles;++t){
        Triangle& t = triangles[t];
        fin.read(reinterpret_cast<char*>tri.normal, 3*sizeof(float));
        fin.read(reinterpret_cast<char*>tri.v0, 3*sizeof(float));
        fin.read(reinterpret_cast<char*>tri.V1, 3*sizeof(float));
        fin.read(reinterpret_cast<char*>tri.v2, 3*sizeof(float));
        char attr[2];
        fin.read(attr,2);
    }
    fin.close();
    return true;
}

static void VoxelReader::voxelizeGrid(Grid3D& grid, const std::vector<Triangle>& triangles){
    std::size_t nx = grid.nx();
    std::size_t ny = grid.ny();
    std::size_t nz = grid.nz();
    //uniform grid assumption todo make it non uniform
    double dx = grid.dx();
    double dy = grid.dx();
    double dz = grid.dx();

    bool isInside = [&](const Vector& voxelCenter, const std::vector<Triangle>& triangles){
        Vector rayDir = {1.0,0,0};
        int crossCount = 0;
        for (const Triangle& tri:triangles){
            if(rayIntersectsTriangle(voxelCenter,rayDir, tri)){
                crossCount++;}
         } 
        return (crossCount%2)==1; 
    }

    for (std::size_t k = 0; k<nz; ++k){
        for (std::size_t j=0; j<ny; ++j){    
            for (std::size_t i=0; i<nx; ++i){
                double x = (i+0.5)*dx;
                double y = (j+0.5)*dy;
                double z = (k+0.5)*dz;
                Vector origin = {x,y,z};
                grid.cellType(i,j,k) = (isInside(origin, triangles))?CellType::INTERIOR:CellType::SOLID;                    
            }
        }
    }
}

static void VoxelReader::voxelizePatch(Grid3D& grid, const std::vector<Triangle>& triangles,const FaceType faceType){
    std::size_t nx = grid.nx();
    std::size_t ny = grid.ny();
    std::size_t nz = grid.nz();
    //uniform grid assumption todo make it non uniform
    double dx = grid.dx();
    double dy = grid.dx();
    double dz = grid.dx();

    for(Triange& tri:triangles){
        float xmin = std::min(tri.v0[0],tri.v1[0],tri.v2[0]);
        float xmax = std::max(tri.v0[0],tri.v1[0],tri.v2[0]);
        float ymin = std::min(tri.v0[1],tri.v1[1],tri.v2[1]);
        float ymax = std::max(tri.v0[1],tri.v1[1],tri.v2[1]);
        float zmin = std::min(tri.v0[2],tri.v1[2],tri.v2[2]);
        float zmax = std::max(tri.v0[2],tri.v1[2],tri.v2[2]);
        std::size_t minI =  floor(xmin/dx);
        std::size_t maxI = min(floor(xmax/dx),nx-1);
        std::size_t minJ =  floor(ymin/dy);
        std::size_t maxJ = max(floor(ymax/dy),ny-1);
        std::size_t minK =  floor(zmin/dz);
        std::size_t maxK = min(floor(zmax/dz),nz-1);
        for (std::size_t k = minK; k<=maxK; ++k){
            for (std::size_t j=minJ; j<maxJ; ++j){    
                for (std::size_t i=minI; i<maxI; ++i){
                    double x = (i+0.5)*dz;
                    double y = (j+0.5)*dy;
                    double z = (k+0.5)*dz;
                    if (grid.cellType(i,j,k) == CellType::BOUNDARY && distanceFromCentroid(x, y, z, tri)<dx)){
                        {
                            grid.faceType(i,j,k) = faceType
                        } 
                    }                 
                }
            }
        }}
}
static double VoxelReader::distanceFromCentroid(double x, double y, double z, const Triangle& tri)
{
    double centX = (tri.v0[0] + tri.v1[0] + tri.v2[0])/3.0;
    double centY = (tri.v0[1] + tri.v1[1] + tri.v2[1])/3.0;
    double centZ = (tri.v0[2] + tri.v1[2] + tri.v2[2])/3.0;
    double distx = centX-x;
    double disty = centy-y;
    double distz = centZ-z;
    return std::sqrt(distx*distx + disty*disty + distz*distz);
}
static bool VoxelReader::rayIntersectsTriangle(const Vector& origin, const Vector& direction, const Triangle& tri) {
    const double epsilon = 1e-8;
    Vector normal = tri.normal;
    double t = normal.dot(tri.v0-origin)/(normal.dot(direction));
    Vector p = origin + direction*t;

    Vector e1 = tri.v1 - tri.v0;
    Vector e2 = tri.v2 - tri.v0;
    Vector vp = p - tri.v0;

    double dot00 = e2.dot(e2);
    double dot01 = e2.dot(e1);
    double dot02 = e2.dot(vp);
    double dot11 = e1.dot(e1);
    double dot12 = e1.dot(vp);

    double bCentDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);

    double u = (dot11 * dot02 - dot01 * dot12) * bCentDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * bCentDenom;

    // inside triangle test
    if (u >= 0.0 && v >= 0.0 && (u + v) <= 1.0)
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
