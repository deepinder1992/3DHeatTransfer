#pragma once
#include "bcType.hpp"
#include <array>
struct SimulationGlobals {
    static constexpr int VERB_LOW    = 1 << 0;
    static constexpr int VERB_MEDIUM = 1 << 1;
    static constexpr int VERB_HIGH   = 1 << 2;
    
    int t = 0;
    int steps = 2;
    int writeInterval = 100000;
    int verbosity = VERB_HIGH;

    double dt = 100;
    double dx = 0.001;
    size_type nx = 8;
    size_type ny = 8;
    size_type nz = 8;
    
    double k = 385; // W/m.K 
    double density = 8960; //kg/m3
    double cp = 385; //J/kg.K
    double alpha = k/(density*cp);

    int maxIters = 100;
    double tol = 1e-6;

    std::array<BCType,6> types = {
                    BCType::Dirichlet, BCType::Neumann, //xmin,max
                    BCType::Neumann, BCType::Dirichlet,    //ymin,max
                    BCType::Dirichlet, BCType::Neumann        //zmin,max
                };  

    std::array<double,6> values = {
                                    100,-10000,
                                    10000,100,
                                    100,-100000};
    // make sure blockdims are power of 2 _best practice
    size_type blockDimX = 8;
    size_type blockDimY = 8;
    size_type blockDimZ = 8;

    

};