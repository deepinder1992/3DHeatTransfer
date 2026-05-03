#pragma once

#include <vector>
#include <array>
#include <random>

#include "simGlobals.hpp"
#include "sparseMatrix.hpp"

// template stays in header
template<typename T>
std::vector<T> randomVector(std::size_t num, T min, T max);

// declarations only (NOT definitions)
extern std::array<BCType,3> types;
extern std::array<double,3> values;

SparseMatrix makeTestMatrix10();

