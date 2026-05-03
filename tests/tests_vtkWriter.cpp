#include <gtest/gtest.h>
#include "../include/grid.hpp"
#include "../include/outputWriter.hpp"

TEST(VTKWriterTest, WritesCorrectValues) {

    std::string dir = "../VTKOutput";
    std::string prefix = "temperature";

    VTKWriter vtkWriter(dir, prefix);;

    Grid3D grid(3,3,3,0.1);

    grid(0,0,0) = 123.0;
    grid(1,0,0) = 456.0;

    vtkWriter.write(grid, 3);

    std::ifstream file(dir + "/" + prefix + "_step_3.vti");

    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());

    EXPECT_NE(content.find("123"), std::string::npos);
    EXPECT_NE(content.find("456"), std::string::npos);
    std::filesystem::remove_all(dir);
}