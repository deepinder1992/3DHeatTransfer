
cmake_minimum_required(VERSION 3.18)

project(Heat3D LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found! Enabling flags.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Default build type = Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set compiler flags per build type
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -Wall -Wextra")

# CUDA debug flags
set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -g -O0")

# Source files
set(SOURCES
    src/main.cpp
    src/grid.cpp
    src/solverCPU.cpp
    src/solverCUDA.cpp
    src/linearAlgebraCPU.cpp
    src/boundaryConditions.cpp
    src/outputWriter.cpp
    src/linearAlgebraGPU.cu
    src/kernel.cu
)

# Create executable
add_executable(heat3d ${SOURCES})

# Include directories
target_include_directories(heat3d PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    /usr/local/cuda/include
)
