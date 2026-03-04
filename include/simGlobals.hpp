#pragma once
#include "bcType.hpp"
#include <array>
struct SimulationGlobals {
    static constexpr int VERB_LOW    = 1 << 0;
    static constexpr int VERB_MEDIUM = 1 << 1;
    static constexpr int VERB_HIGH   = 1 << 2;
    
    int t = 0;
    int steps = 1;
    int writeInterval = 100000;
    double globalTol = 1e-8;
    int verbosity = VERB_HIGH;

    double dt = 1;
    double lx = 1; // 10 cm =ly,lz
   // double ly = 0.1; // 10 cm
    //double lz = 0.1; // 10 cm

 
    size_type nx = 60;
    size_type ny = nx;
    size_type nz = nx;

    double dx = lx/nx; // = dy,dz
   // double dy = ly/ny;
   // double dz = lz/nz;


    
    double k = 385; // W/m.K 
    double density = 8960; //kg/m3
    double cp = 385; //J/kg.K
    double alpha = k/(density*cp);

    int maxIters = 50;
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
    size_type blockDimX = 8;
    size_type blockDimY = 8;
    size_type blockDimZ = 8;
    //size_type blockDim1D = blockDimX*blockDimY*blockDimZ;
    //int numSMs = 16; // can vary based on GPU model
    //size_type gridDim1D = 6*numSMs;
    //debugging
    mutable int totalIters = 0;
    

};