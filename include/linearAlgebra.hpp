#pragma once
#include <cmath>
#include "sparseMatrix.hpp"
#include <iostream>

double dot (const std::vector<double>& a, const std::vector<double>& b){
    double sum = 0.0;
    for (size_type i = 0; i <a.size();++i)
            sum+=a[i]*b[i];
    return sum;
}


void SparseMultiply (const SparseMatrix& A, const double* x, double* y){
    
    const auto& values = A.values();
    const auto& cols   = A.colIndex();
    const auto& rowPtr = A.rowPtr();

    for (size_type row = 0; row < A.rows(); ++row){
        double sum = 0.0;       
        for (size_type idx = rowPtr[row]; idx <rowPtr[row+1]; ++idx){
            sum+=values[idx]*x[cols[idx]];
        }
        y[row] = sum;
    }
}

void conjugateGradient(const SparseMatrix& A,
                       const std::vector<double>& b,
                       std::vector<double>& x,
                         const SimulationGlobals& globs)
{
    int N = b.size();
    std::vector<double> r(N), p(N), Ap(N);

    SparseMultiply(A, x.data(), Ap.data());
    for (int i = 0; i < N; ++i)
        r[i] = b[i] - Ap[i];

    p = r;
    double rsold = dot(r, r);

    for (int iter = 0; iter < globs.maxIters; ++iter)
    {
        SparseMultiply(A, p.data(), Ap.data());
        double alpha = rsold / dot(p, Ap);

        for (int i = 0; i < N; ++i){
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];}

        double rsnew = dot(r, r);
        double sqrtRsnew = std::sqrt(rsnew);
        if (globs.verbosity & SimulationGlobals::VERB_HIGH){
            std::cout << "     Step:: "<<globs.t+1<<" Iter:  "<< iter<< "  Err:  "<<sqrtRsnew << std::endl;
            ++globs.totalIters;}
        if (sqrtRsnew < globs.tol) break;

        for (int i = 0; i < N; ++i)
            p[i] = r[i] + (rsnew / rsold) * p[i];

        rsold = rsnew;
    }
}
