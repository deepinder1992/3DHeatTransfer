#include "grid.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iostream>

class OutputWriter{
    public:
        virtual ~OutputWriter() =default;
        virtual void write(const Grid3D& grid, const int& t)= 0;
};


class VTKWriter: public OutputWriter{
    public:
    VTKWriter(const std::string directory,const std::string& prefix):directory_(directory),prefix_(prefix)
             {std::filesystem::create_directories(directory_);}
    void write(const Grid3D& grid, const int& t) override;

    private:
        std::string directory_,prefix_;
};
