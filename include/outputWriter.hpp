#include "grid.hpp"
#include <string>
#include <fstream>
class OutputWriter{
    public:
        virtual ~OutputWriter() =default;
        virtual void write(const Grid3D& grid, int step)= 0;
};

class BinaryWriter: OutputWriter{
    public:
        BinaryWriter(const std::string& prefix):prefix_(prefix){}
        
        void write(const Grid3D& grid, int step) override;
    
        private:
            std::string prefix_;
};