#pragma once
#include<vector>
#include<cassert>
#include<cstddef>
#include "cuda_runtime.h"
#include "simGlobals.hpp"

using size_type = std::size_t;
class Grid3D{
    public:
        Grid3D(size_type nx, size_type ny, size_type nz);
        
        double& operator()(size_type i, size_type j, size_type k);
        const double&  operator()(size_type i, size_type j, size_type k) const;

        CellType& cellType(size_type i, size_type j,size_type k);
        const CellType& cellType(size_type i, size_type j,size_type k) const;


        size_type nx() const noexcept {return nx_;}
        size_type ny() const noexcept {return ny_;}
        size_type nz() const noexcept {return nz_;}

        size_type size() const noexcept{return data_.size();}
        double* data() noexcept {return data_.data();} 
        const double* data() const noexcept{return data_.data();}
        
        void fill (double value);
        
        void detectBoundaries();
    private:
        size_type index(size_type i, size_type j, size_type k) const noexcept;
        
        size_type nx_, ny_,nz_;

        std::vector<double> data_;

        std::vector<CellType> cellType_;

};