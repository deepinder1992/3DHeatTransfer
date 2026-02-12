#include "outputWriter.hpp"

void BinaryWriter::write(const Grid3D& grid, int step){
    std::ofstream out(directory_+"/"+ prefix_ +"_"+ std::to_string(step) + ".bin", std::ios::binary);
    out.write(reinterpret_cast<const char*>(grid.data()), grid.size()*sizeof(double));
}

void VTKWriter::write(const Grid3D& grid, int step){
    std::ofstream file;
    std::stringstream filename;
    filename <<directory_<<"/"<<prefix_<<"_step_"<<step <<".vti";
    file.open(filename.str());
    
    if (!file.is_open()){
        std::cerr<<"Failed to Open file for writing: "<< filename.str()<<std::endl;
    }

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <ImageData WholeExtent=\"0 " << grid.nx()-1 << " 0 " << grid.ny()-1 << " 0 " << grid.nz()-1<<"\"" <<
                                                    " Origin="<<"\"0 0 0\""<< 
                                                    " Spacing="<<"\"1 1 1\""<< ">\n";
    file << "    <Piece Extent=\"0 " << grid.nx()-1 << " 0 " << grid.ny()-1 << " 0 " << grid.nz()-1 << "\">\n";
    file << "      <PointData Scalars=\"temperature\">\n";

    file << "        <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">\n";
    for (size_t k = 0; k < grid.nz(); ++k) {
        for (size_t j = 0; j < grid.ny(); ++j) {
            for (size_t i = 0; i < grid.nx(); ++i) {
                file << grid(i, j, k) << " ";
            }
        }
    }
    file << "\n        </DataArray>\n";
    file << "      </PointData>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "</VTKFile>\n";

    file.close();
    }

                    

