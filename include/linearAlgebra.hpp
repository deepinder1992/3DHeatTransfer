#pragma once
#include <cmath>
#include "sparseMatrix.hpp"
#include <iostream>
#include "simGlobals.hpp"

class LinearAlgebra{
    public:         
         LinearAlgebra(int maxIters):_maxIters(maxIters){};

        double dot (const std::vector<double>& a, const std::vector<double>& b);

        void SparseMultiply (const SparseMatrix& A, const double* x, double* y);

        void conjugateGradient(const SparseMatrix& A, const std::vector<double>& b,
                                std::vector<double>& x, const SimulationGlobals& globs);

        void implicitJacobiCPU(std::size_t nx, std::size_t ny, std::size_t nz, const double coeff_, double& maxerror,
                                 Grid3D* oldGrid, Grid3D* newGrid, const Grid3D& current);

        int maxIters() const noexcept{return _maxIters;}

        //in case of non convergnce increase max iters but cap at 2000 to avoid infinte loop
        void adjustMaxItersIfNeeded(int iter){
                    if (iter==(_maxIters-1) && _maxIters<2000){
                        _maxIters = static_cast<int> (_maxIters*1.5);
                        std::cout<< "Inner Solution did not converge increasing maxIters to " << _maxIters<<std::endl;
                    } }
                        
private:
    int _maxIters = 0;
};