
# Minimum required version
cmake_minimum_required(VERSION 3.18)

# Project name and languages
project(Heat3D LANGUAGES CXX CUDA)

# -----------------------------
# C++ Standard
# -----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------
# Default build type = Debug
# -----------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -O0")
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -lineinfo")
endif()


# -----------------------------
# Include directories
# -----------------------------
include_directories(${PROJECT_SOURCE_DIR}/include)


set(SOURCES
    src/main.cpp
    src/grid.cpp
    src/solverCPU.cpp
    src/boundaryConditions.cpp
)

add_executable(heat3d ${SOURCES})