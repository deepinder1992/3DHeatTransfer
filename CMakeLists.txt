cmake_minimum_required(VERSION 3.18)

project(Heat3D LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect packages
find_package(CUDAToolkit REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "CUDA Toolkit found! Enabling flags.")
endif()

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
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")

# CUDA debug flags
set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -g -O0")

set(SHARED_SOURCES
    src/grid.cpp
    src/boundaryConditions.cpp
    src/solverCPU.cpp
    src/linearAlgebraCPU.cpp
    src/kernel.cu
    src/voxelReader.cpp
)

set_source_files_properties(
    src/boundaryConditions.cpp
    PROPERTIES LANGUAGE CUDA
)
# Source files
set(SOURCES
src/main.cpp
    src/solverCUDAStencil.cpp
    src/solverCUDAMatrix.cpp
    src/outputWriter.cpp
    src/linearAlgebraGPU.cu
    ${SHARED_SOURCES}  
)

# Create executable
add_executable(heat3d ${SOURCES})

# Include directories
target_include_directories(heat3d PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    /usr/local/cuda/include
)

### tests
add_executable(test_heat
    tests/tests_main.cu
	${SHARED_SOURCES}  
	)
target_include_directories(test_heat PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(test_heat PRIVATE
    OpenMP::OpenMP_CXX
    CUDA::cudart
)
