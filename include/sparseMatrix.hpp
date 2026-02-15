#pragma once
#include <vector>
#include "grid.hpp"

class SparseMatrix {
    public:
        SparseMatrix(size_type n): _nRows(n),_rowPtr(n+1,0){}

        size_type rows()const {return _nRows;}

        std::vector<double>& values() {return _values;}
        std::vector<size_type>& colIndex() {return _colIndex;}
        std::vector<size_type>& rowPtr() {return _rowPtr;}

        const std::vector<double>& values() const {return _values;}
        const std::vector<size_type>& colIndex() const{return _colIndex;}
        const std::vector<size_type>& rowPtr() const {return _rowPtr;}
    
    private:
        size_type _nRows;
        std::vector<double> _values;
        std::vector<size_type> _colIndex;
        std::vector<size_type> _rowPtr;
};