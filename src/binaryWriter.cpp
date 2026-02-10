#include "outputWriter.hpp"

void BinaryWriter::write(const Grid3D& grid, int step){
    std::ofstream out(prefix_ +"_"+ std::to_string(step) + ".bin", std::ios::binary);
    out.write(reinterpret_cast<const char*>(grid.data()), grid.size()*sizeof(double));
}
