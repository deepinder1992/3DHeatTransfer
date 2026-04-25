cmake_minimum_required(VERSION 3.20)

option(ENABLE_CUDA "Enable CUDA support for GPU solvers" OFF)

project(HeatTransfer3D LANGUAGES CXX)

if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_definitions(-DENABLE_CUDA)
    message(STATUS "CUDA enabled")
else()
    message(STATUS "CUDA disabled - building CPU-only version")
endif()

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Find dependencies
find_package(OpenMP REQUIRED)

# Source files
set(COMMON_SOURCES
    src/grid.cpp
    src/boundaryConditions.cpp
    src/solverCPU.cpp
    src/linearAlgebraCPU.cpp
    src/voxelReader.cpp
    src/outputWriter.cpp
)

set(SOURCES
    src/main.cpp
    ${COMMON_SOURCES}
)

if(ENABLE_CUDA)
    list(APPEND SOURCES
        src/solverCUDAStencil.cpp
        src/solverCUDAMatrix.cpp
        src/linearAlgebraGPU.cu
        src/kernel.cu
    )
endif()

# Main executable
add_executable(heat3d ${SOURCES})

target_include_directories(heat3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(heat3d PRIVATE OpenMP::OpenMP_CXX)

if(ENABLE_CUDA)
    target_link_libraries(heat3d PRIVATE CUDA::cudart)
    set_target_properties(heat3d PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

# Compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(heat3d PRIVATE -Wall -Wextra -O0)
else()
    target_compile_options(heat3d PRIVATE -O3 -march=native -Wall -Wextra)
endif()

# Tests
add_executable(test_heat
    tests/tests_main.cpp
    ${COMMON_SOURCES}
)

target_include_directories(test_heat PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_heat PRIVATE OpenMP::OpenMP_CXX)

if(ENABLE_CUDA)
    target_sources(test_heat PRIVATE tests/tests_main.cu)
    target_link_libraries(test_heat PRIVATE CUDA::cudart)
endif()

enable_testing()
add_test(NAME heat_tests COMMAND test_heat)
