#pragma once
#include "bcType.hpp"
#include <array>
struct SimulationGlobals {
    static constexpr int VERB_LOW    = 1 << 0;
    static constexpr int VERB_MEDIUM = 1 << 1;
    static constexpr int VERB_HIGH   = 1 << 2;
    
    int t = 0;
    int steps = 10000;
    int writeInterval = 100;
    int verbosity = VERB_LOW;

    double dt = 1;
    double dx = 1;
    
    double k = 1;//385; // W/m.K 
    double density =1;// 8960; //kg/m3
    double cp = 1;//385; //J/kg.K
    double alpha = k/(density*cp);

    int maxIters = 50;
    double tol = 1e-6;

    std::array<BCType,6> types = {
                    BCType::Dirichlet, BCType::Neumann, //xmin,max
                    BCType::Neumann, BCType::Dirichlet,    //ymin,max
                    BCType::Dirichlet, BCType::Neumann        //zmin,max
                };  

    std::array<double,6> values = {
                                    100,-10,
                                    0,100,
                                    100,0};
    

};