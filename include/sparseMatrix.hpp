#pragma once
#include <vector>
#include "grid.hpp"

class SparseMatrix {
    public:
        SparseMatrix(std::size_t n): _nRows(n),_rowPtr(n+1,0){}

        std::size_t rows()const {return _nRows;}

        std::vector<double>& values() {return _values;}
        std::vector<std::size_t>& colIndex() {return _colIndex;}
        std::vector<std::size_t>& rowPtr() {return _rowPtr;}

        const std::vector<double>& values() const {return _values;}
        const std::vector<std::size_t>& colIndex() const{return _colIndex;}
        const std::vector<std::size_t>& rowPtr() const {return _rowPtr;}
    
    private:
        std::size_t _nRows;
        std::vector<double> _values;
        std::vector<std::size_t> _colIndex;
        std::vector<std::size_t> _rowPtr;
};