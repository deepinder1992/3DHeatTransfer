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


if(ENABLE_CUDA)
    list(APPEND COMMON_SOURCES
        src/cudaSrc/solverCUDAStencil.cu
        src/cudaSrc/solverCUDAMatrix.cu
        src/cudaSrc/linearAlgebraGPU.cu
        src/cudaSrc/kernel.cu
        src/cudaSrc/boundaryConditions.cu
    )
endif()

# Main executable
add_executable(heat3d 
    src/main.cpp
    ${COMMON_SOURCES})

target_include_directories(heat3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(heat3d PRIVATE OpenMP::OpenMP_CXX)

if(ENABLE_CUDA)
    target_link_libraries(heat3d PRIVATE CUDA::cudart)
    set_target_properties(heat3d PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
endif()

# Compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(heat3d PRIVATE -Wall -Wextra -O0)
else()
    target_compile_options(heat3d PRIVATE -O3 -march=native -Wall -Wextra)
endif()

#tests
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings  
# https://google.github.io/googletest/quickstart-cmake.html
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Tests
add_executable(test_heat
    tests/tests_main.cpp
    tests/tests_grid.cpp
    tests/tests_VoxelReader.cpp
    tests/tests_cpuLinearAlgebra.cpp
    tests/tests_sparseMatrix.cpp
    tests/tests_boundaryConditions.cpp
    tests/test_stencilFulltest.cpp
    tests/test_matrixFulltest.cpp
    tests/tests_vtkWriter.cpp
    tests/tests_utils.cpp
    tests/tests_analytical.cpp
    ${COMMON_SOURCES}
)

target_include_directories(test_heat PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_heat PRIVATE OpenMP::OpenMP_CXX
                        gtest
                        gtest_main)

if(ENABLE_CUDA)
    target_sources(test_heat  
    PRIVATE tests/tests_cudaLinAlgebra.cu)
    target_link_libraries(test_heat PRIVATE CUDA::cudart)
endif()

enable_testing()
include(GoogleTest)
gtest_discover_tests(test_heat)
