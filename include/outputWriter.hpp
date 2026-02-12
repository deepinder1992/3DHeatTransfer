#include "grid.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iostream>

class OutputWriter{
    public:
        virtual ~OutputWriter() =default;
        virtual void write(const Grid3D& grid, int step)= 0;
};

class BinaryWriter: public OutputWriter{
    public:
        BinaryWriter(const std::string directory,const std::string& prefix):directory_(directory),prefix_(prefix)
                    {std::filesystem::create_directories(directory_);}
        
        void write(const Grid3D& grid, int step) override;
    
    private:
        std::string prefix_,directory_;
};

class VTKWriter: public OutputWriter{
    public:
    VTKWriter(const std::string directory,const std::string& prefix):directory_(directory),prefix_(prefix)
             {std::filesystem::create_directories(directory_);}
    void write(const Grid3D& grid, int step) override;

    private:
        std::string prefix_,directory_;
};
