#pragma once
#include <array>
#include <string>

enum class SolverType
{
    CPU_STENCIL = 1,
    CPU_MATRIX  = 2,
    CUDA_STENCIL = 3,
    CUDA_MATRIX  = 4
};

enum class BCType{
    Dirichlet =1,
    Neumann = 2
};

enum class CellType
{
    SOLID = 0,
    INTERIOR = 1,
    BOUNDARY = 2
};
enum class FaceType
{
    INLET = 0,
    OUTLET = 1,
    WALL = 2,
    NONE = 4
};
struct SimulationGlobals {
    static constexpr int VERB_LOW    = 1 << 0;
    static constexpr int VERB_MEDIUM = 1 << 1;
    static constexpr int VERB_HIGH   = 1 << 2;

    SolverType solver = SolverType::CPU_STENCIL; // default solver
    std::string stlFileloc  = "../stlFiles/cylinder.stl"; //todo make it os agnostic
    
    int t = 0;
    int steps = 500;
    int writeInterval = 10;
    double globalTol = 1e-8;
    int verbosity = VERB_LOW;

    double dt = 50;
    //double lx = 60; //  =ly,lz

    std::size_t nx = 50;
    std::size_t ny = nx;
    std::size_t nz = nx;

    double dx = 0.1; // = dy,dz
   // double dy = ly/ny;
   // double dz = lz/nz;


    
    double k = 385; // W/m.K 
    double density = 8960; //kg/m3
    double cp = 385; //J/kg.K
    double alpha = k/(density*cp);

    mutable int maxIters = 50;
    double tol = 1e-6;

    std::array<BCType,3> types = {
                    BCType::Neumann, //inlet
                    BCType::Dirichlet,    //outlet
                    BCType::Neumann        //wall
                };  

    std::array<double,3> values = { 50000, //inlet
                                    100, //outlet
                                    -50000}; //wall
    // make sure blockdims are power of 2 _best practice
    std::size_t blockDimX = 8;
    std::size_t blockDimY = 8;
    std::size_t blockDimZ = 8;
    mutable int totalIters = 0;   

};

