#pragma once

#include<vector>
#include<cassert>
#include<cstddef>

class Grid3D{
    public:
        using value_type = double;
        using size_type = std::size_t;


        Grid3D(size_type nx, size_type ny, size_type nz);

        value_type& operator()(size_type i, size_type j, size_type k);
        const value_type&  operator()(size_type i, size_type j, size_type k) const;

        size_type nx() const noexcept {return nx;}
        size_type ny() const noexcept {return ny;}
        size_type nz() const noexcept {return nz;}

        size_type size() const noexcept{return data_.Size();}
        value_type* data() noexcept {return data_.data();} 
        const value_type* data() const noexcept{return data_.data();}

        void fill (value_type value);
    
    private:
        size_type index(size_type i, size_type j, size_type k) noexcept;
        
        size_type nx_, ny_,nz_;

        std::vector<value_type> data_;

}