#include "linearAlgebra.hpp"

double LinearAlgebra::dot (const std::vector<double>& a, const std::vector<double>& b){
    double sum = 0.0;
    for (std::size_t i = 0; i <a.size();++i)
            sum+=a[i]*b[i];
    return sum;
}


void LinearAlgebra::SparseMultiply (const SparseMatrix& A, const double* x, double* y){
    
    const auto& values = A.values();
    const auto& cols   = A.colIndex();
    const auto& rowPtr = A.rowPtr();

    for (std::size_t row = 0; row < A.rows(); ++row){
        double sum = 0.0;       
        for (std::size_t idx = rowPtr[row]; idx <rowPtr[row+1]; ++idx){
            sum+=values[idx]*x[cols[idx]];
        }
        y[row] = sum;
    }
}

void LinearAlgebra::conjugateGradient(const SparseMatrix& A,
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
    for (int iter = 0; iter < _maxIters; ++iter)
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
        adjustMaxItersIfNeeded(iter);
    }
}

void LinearAlgebra::implicitJacobiCPU(std::size_t nx, std::size_t ny, std::size_t nz, const double coeff_, double& maxErr,
                                     Grid3D* oldGrid, Grid3D* newGrid, const Grid3D& current){
        #pragma omp parallel for collapse(3) reduction(max:maxErr)
        for (std::size_t k = 1; k < nz-1; ++k){
            for (std::size_t j = 1; j < ny-1; ++j ){
                for (std::size_t i = 1; i < nx-1 ; ++i){
                    if(current.cellType(i,j,k)!=CellType::INTERIOR) continue;
                    const double rhs  = current(i,j,k);
                    
                    const auto& old = *oldGrid;

                    const double sum = old(i+1,j,k)+ old(i-1,j,k)
                                                + old(i,j+1,k)+ old(i,j-1,k)
                                                    + old(i,j,k+1)+ old(i,j,k-1);
                    
                    const double newVal = (rhs + coeff_*sum)/(1+6*coeff_);
                    
                    maxErr = std::max(maxErr, std::abs(newVal-old(i,j,k)));
                    (*newGrid)(i,j,k) = newVal;
                }
            }        
        }
}
    
