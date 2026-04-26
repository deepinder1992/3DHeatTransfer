#pragma once

#include <vector>
#include <random>

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