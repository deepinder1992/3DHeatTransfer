#pragma once

#include<vector>
#include<cassert>
#include<cstddef>
using size_type = std::size_t;
class Grid3D{
    public:
        using value_type = double;

        Grid3D(size_type nx, size_type ny, size_type nz);

        value_type& operator()(size_type i, size_type j, size_type k);
        const value_type&  operator()(size_type i, size_type j, size_type k) const;

        size_type nx() const noexcept {return nx_;}
        size_type ny() const noexcept {return ny_;}
        size_type nz() const noexcept {return nz_;}

        size_type size() const noexcept{return data_.size();}
        value_type* data() noexcept {return data_.data();} 
        const value_type* data() const noexcept{return data_.data();}

        void fill (value_type value);
    
    private:
        size_type index(size_type i, size_type j, size_type k) const noexcept;
        
        size_type nx_, ny_,nz_;

        std::vector<value_type> data_;

};