#include "./include/tests_utils.hpp"

// template implementation
template<typename T>
std::vector<T> randomVector(std::size_t num, T min, T max) {
    std::vector<T> randVec(num);
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<T> dist(min, max);

    for (auto &x : randVec) {
        x = dist(gen);
    }

    return randVec;
}

// explicit instantiations (important!)
template std::vector<double> randomVector<double>(std::size_t, double, double);
template std::vector<float> randomVector<float>(std::size_t, float, float);

// matrix function
SparseMatrix makeTestMatrix10() {
    SparseMatrix mat(10);

    mat.rowPtr() = {0,2,5,8,11,14,17,20,23,26,29};

    mat.values() = {
        2,1,
        -1,2,1,
        -1,2,1,
        -1,2,1,
        -1,2,1,
        -1,2,1,
        -1,2,1,
        -1,2,1,
        -1,2,1,
        -1,2
    };

    mat.colIndex() = {
        0,1,
        0,1,2,
        1,2,3,
        2,3,4,
        3,4,5,
        4,5,6,
        5,6,7,
        6,7,8,
        7,8,9,
        8,9
    };

    return mat;
}

