#pragma once
#include <array>

enum class SolverType
{
    CPU_STENCIL = 1,
    CPU_MATRIX  = 2,
    CUDA_STENCIL = 3,
    CUDA_MATRIX  = 4
};

enum class BCType{
    Dirichlet,
    Neumann
};

struct SimulationGlobals {
    static constexpr int VERB_LOW    = 1 << 0;
    static constexpr int VERB_MEDIUM = 1 << 1;
    static constexpr int VERB_HIGH   = 1 << 2;

    SolverType solver = SolverType::CPU_MATRIX; // default solver
    
    int t = 0;
    int steps = 1000;
    int writeInterval = 100000;
    double globalTol = 1e-8;
    int verbosity = VERB_LOW;

    double dt = 50;
    double lx = 1; // 10 cm =ly,lz


 
    std::size_t nx = 30;
    std::size_t ny = nx;
    std::size_t nz = nx;

    double dx = lx/nx; // = dy,dz
   // double dy = ly/ny;
   // double dz = lz/nz;


    
    double k = 385; // W/m.K 
    double density = 8960; //kg/m3
    double cp = 385; //J/kg.K
    double alpha = k/(density*cp);

    mutable int maxIters = 50;
    double tol = 1e-6;

    std::array<BCType,6> types = {
                    BCType::Dirichlet, BCType::Neumann, //xmin,max
                    BCType::Neumann, BCType::Dirichlet,    //ymin,max
                    BCType::Dirichlet, BCType::Neumann        //zmin,max
                };  

    std::array<double,6> values = {
                                    100,-50000,
                                    50000,100,
                                    100,-50000};
    // make sure blockdims are power of 2 _best practice
    std::size_t blockDimX = 8;
    std::size_t blockDimY = 8;
    std::size_t blockDimZ = 8;
    //std::size_t blockDim1D = blockDimX*blockDimY*blockDimZ;
    //int numSMs = 16; // can vary based on GPU model
    //std::size_t gridDim1D = 6*numSMs;
    //debugging
    mutable int totalIters = 0;
    

};

