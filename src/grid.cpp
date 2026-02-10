#include "grid.hpp"


Grid3D::Grid3D(size_type nx, size_type ny, size_type nz)
:nx_(nx),ny_(ny),nz_(nz), data_(nx*ny*nz)
    {
        assert(nx>0 && ny > 0 && nz >0);
    }

Grid3D::value_type& Grid3D::operator()(size_type i, size_type j, size_type k)
    {  assert(i<nx_ && j< ny_ && k <nz_);
        return data_[index(i,j,k)];
    }


const Grid3D::value_type& Grid3D::operator()(size_type i, size_type j, size_type k) const
    {  assert(i<nx_ && j< ny_ && k <nz_);
        return data_[index(i,j,k)];
    }


void Grid3D::fill(value_type value)
    {
        std::fill(data_.begin(),data_.end(),value);
    }

size_type Grid3D::index(size_type i, size_type j, size_type k) const noexcept
    {
        return i+ nx_*(j+ ny_*k);
    }