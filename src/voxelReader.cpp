#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "voxelReader.hpp"

bool rayIntersectsTriangle(const Vector& origin, const Vector& direction, const Triangle& tri) {
    constexpr double epsilon = 1e-8;
    constexpr double epsilon2 = 1e-7;
    Vector normal = tri.normal;
    double nDenom = normal.dot(direction);
    if(std::abs(nDenom) < epsilon){return false;}// with in parallel threshold
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
    if (std::abs(bCentDenom) <= epsilon){return false;} //degenrate triangle

    double u = (dot11 * dot02 - dot01 * dot12) / bCentDenom;
    double v = (dot00 * dot12 - dot01 * dot02) / bCentDenom;

    // inside triangle test
    if (u >= -epsilon && v >= -epsilon && (u + v) <= 1.0+epsilon)
        return true;

    return false;
}

float VoxelReader::bBoxMinX = 0.0;
float VoxelReader::bBoxMaxX = 0.0;
float VoxelReader::bBoxMinY = 0.0;
float VoxelReader::bBoxMaxY = 0.0;
float VoxelReader::bBoxMinZ = 0.0;
float VoxelReader::bBoxMaxZ = 0.0;

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
        static_assert(sizeof(Vector) == 3 * sizeof(float), "Vector must have no padding");
        std::cout << "Reading Stl file..\n";
        loadBinaryStl(fileName, triangles);

        std::cout <<"Started Voxelizing..\n";
        grid.adjustGrid(maxLen());
        shiftTriangles(triangles,grid.gridCent());

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
        shiftTriangles(inletTriangles,grid.gridCent());

        voxelizePatch(grid, inletTriangles, FaceType::INLET);
    

        std::cout << "Applying outlet patch..\n";
        std::vector<Triangle> outletTriangles;
        loadBinaryStl(outletFile, outletTriangles);
        shiftTriangles(outletTriangles,grid.gridCent());

        voxelizePatch(grid, outletTriangles, FaceType::OUTLET);
        
        std::cout << "Applying wall patch..\n";
        std::vector<Triangle> wallTriangles;
        loadBinaryStl(wallFile, wallTriangles);
        shiftTriangles(wallTriangles,grid.gridCent());
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
    
    fin.read(reinterpret_cast<char*>(&numTriangles), sizeof(uint32_t));
    triangles.resize(numTriangles);
    for (uint32_t t=0; t<numTriangles;++t){
        Triangle& tri = triangles[t];
        fin.read(reinterpret_cast<char*>(&tri.normal), 3*sizeof(float));
        fin.read(reinterpret_cast<char*>(&tri.v0), 3*sizeof(float));
        fin.read(reinterpret_cast<char*>(&tri.v1), 3*sizeof(float));
        fin.read(reinterpret_cast<char*>(&tri.v2), 3*sizeof(float));
        char attr[2];
        fin.read(attr,2);
        boundingBox(tri, bBoxMinX, bBoxMaxX, bBoxMinY, bBoxMaxY, bBoxMinZ, bBoxMaxZ);
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
    int kkkk = 0;
    auto isInside = [&](const Vector& voxelCenter, const std::vector<Triangle>& triangles){
        Vector rayDir = {1.0,0,0};
        int crossCount = 0;
        for (const Triangle& tri:triangles){
            if(rayIntersectsTriangle(voxelCenter,rayDir, tri)){
                crossCount++;}
         } 
        if((crossCount%2)==1)kkkk+=1;
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
    std::cout<<kkkk<<std::endl;
}

void VoxelReader::voxelizePatch(Grid3D& grid, const std::vector<Triangle>& triangles,const FaceType faceType){
    std::size_t nx = grid.nx();
    std::size_t ny = grid.ny();
    std::size_t nz = grid.nz();
    //uniform grid assumption todo make it non uniform
    double dx = grid.dx();
    double dy = grid.dx();
    double dz = grid.dx();
    int kkkk = 0;

    for(const Triangle& tri:triangles){
        float xmin = 0.0f, xmax = 0.0f, ymin = 0.0f, ymax = 0.0f, zmin = 0.0f, zmax = 0.0f;
        boundingBox(tri, xmin, xmax, ymin, ymax, zmin, zmax); 
        std::size_t minI = std::max(static_cast<std::size_t>(floor(xmin / dx)), size_t(0));
        std::size_t maxI = std::min(static_cast<std::size_t>(ceil(xmax / dx)), nx - 1);
        std::size_t minJ = std::max(static_cast<std::size_t>(floor(ymin / dy)), size_t(0));
        std::size_t maxJ = std::min(static_cast<std::size_t>(ceil(ymax / dy)), ny - 1);
        std::size_t minK = std::max(static_cast<std::size_t>(floor(zmin / dz)), size_t(0));
        std::size_t maxK = std::min(static_cast<std::size_t>(ceil(zmax / dz)), nz - 1);
        for (std::size_t k = minK; k<=maxK; ++k){
            for (std::size_t j=minJ; j<=maxJ; ++j){    
                for (std::size_t i=minI; i<=maxI; ++i){
                    float x = (i+0.5)*dx;
                    float y = (j+0.5)*dy;
                    float z = (k+0.5)*dz;
                    if (grid.cellType(i,j,k) == CellType::BOUNDARY //&& grid.faceType(i,j,k) == FaceType::NONE
                         && isInterSecting(x, y, z, dx/2.0,tri )){
                        // && distanceFromCentroid(x, y, z, tri)<=(dx+1e-3f)){
                            grid.faceType(i,j,k) = faceType;
                            grid.cellFaceNormal(i,j,k) = tri.normal;
                            ++kkkk;
                        
                    }                 
                }
            }
        }}
    std::cout<<"iiii"<<kkkk<<std::endl;
}

void VoxelReader::boundingBox(const Triangle& tri, float& minX,float& maxX, float& minY,
                                float& maxY,float& minZ, float& maxZ ){

    const Vector verts[3] = {tri.v0, tri.v1, tri.v2};

    for (int i = 0; i < 3; i++) {

        minX = std::min(minX, verts[i].x);
        minY = std::min(minY, verts[i].y);
        minZ = std::min(minZ, verts[i].z);

        maxX = std::max(maxX, verts[i].x);
        maxY = std::max(maxY, verts[i].y);
        maxZ = std::max(maxZ, verts[i].z);
    }
}

void VoxelReader::shiftTriangles(std::vector<Triangle>& triangles, Vector gridCent) {
    Vector triangleCent = {0.5f*(bBoxMaxX+bBoxMinX), 0.5f*(bBoxMaxY+bBoxMinY), 0.5f*(bBoxMaxZ+bBoxMinZ)};
    Vector shiftVect = gridCent-triangleCent;

    for (Triangle& tri : triangles) {
        tri.v0.x += shiftVect.x;
        tri.v0.y += shiftVect.y;
        tri.v0.z += shiftVect.z;

        tri.v1.x += shiftVect.x;
        tri.v1.y += shiftVect.y;
        tri.v1.z += shiftVect.z;

        tri.v2.x += shiftVect.x;
        tri.v2.y += shiftVect.y;
        tri.v2.z += shiftVect.z;
    }
}

double VoxelReader::maxLen() {
    double lenx = fabs(bBoxMinX-bBoxMaxX);
    double leny = fabs(bBoxMinY-bBoxMaxY);
    double lenz = fabs(bBoxMinZ-bBoxMaxZ);
    
    return std::max(lenx,std::max(leny,lenz));
}

bool VoxelReader::isInterSecting(float boxCentX, float boxCentY, float boxCentZ, float halfSize, const Triangle& tri ){
    //Separating Axis Theorem
    const float eps = 1.0e-3f * halfSize; 
    Vector boxCent = {boxCentX, boxCentY, boxCentZ};

    Vector v0 = tri.v0-boxCent;
    Vector v1 = tri.v1-boxCent;
    Vector v2 = tri.v2-boxCent;

    auto maxSpan = [](float a, float b, float c){return std::max(a,std::max(b,c));};
    auto minSpan = [](float a, float b, float c){return std::min(a,std::min(b,c));};

    if(minSpan(v0.x, v1.x, v2.x)>halfSize+eps) {return false;}
    if(maxSpan(v0.x, v1.x, v2.x)<-halfSize-eps) {return false;}

    if(minSpan(v0.y, v1.y, v2.y)>halfSize+eps) {return false;}
    if(maxSpan(v0.y, v1.y, v2.y)<-halfSize-eps) {return false;}

    if(minSpan(v0.z, v1.z, v2.z)>halfSize+eps) {return false;}
    if(maxSpan(v0.z, v1.z, v2.z)<-halfSize-eps) {return false;}

    Vector e0 = v0-v1;
    Vector e1 = v1-v2;
    Vector e2 = v2-v0;
    //plane vs triangle
    Vector triNorm = e0^e1;
    float planeConst = -triNorm.dot(v0);
    Vector vmin,vmax;
    for(int i =0; i<3; ++i){
        float normComp = (&triNorm.x)[i];
        if(normComp>0.0){
            (&vmin.x)[i] = -halfSize;
            (&vmax.x)[i] = halfSize;

        }
        else{
            (&vmax.x)[i] = -halfSize;
            (&vmin.x)[i] = halfSize;
        }
    }

    if(triNorm.dot(vmin)+planeConst>eps){ return false;}
    if(triNorm.dot(vmax)+planeConst<-eps){ return false;}

    Vector boxAxes[3] = { {1,0,0}, {0,1,0}, {0,0,1} };
    Vector edges[3] = { e0, e1, e2 };

    for (int i = 0; i < 3; ++i) {       
        for (int j = 0; j < 3; ++j) {   
            Vector axis = edges[i] ^ boxAxes[j];  
            if (axis.x == 0.0f && axis.y == 0.0f && axis.z == 0.0f) continue; 

            float triProjMin = std::min({v0.dot(axis), v1.dot(axis), v2.dot(axis)});
            float triProjMax = std::max({v0.dot(axis), v1.dot(axis), v2.dot(axis)});
   
            float r = halfSize * (std::fabs(axis.x) + std::fabs(axis.y) + std::fabs(axis.z));
    
            if (triProjMin > r || triProjMax < -r) return false;
        }
    }
    return true;
}
