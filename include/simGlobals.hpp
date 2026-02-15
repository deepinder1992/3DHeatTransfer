#pragma once

struct SimulationGlobals {
    int t = 0;
    int steps = 100;
    int writeInterval = 1;
    int verbosity = 1;

    double dt = 1.0;
    double alpha = 0.1;
    double dx = 0.1;

    static constexpr int VERB_LOW    = 1 << 0;
    static constexpr int VERB_MEDIUM = 1 << 1;
    static constexpr int VERB_HIGH   = 1 << 2;
};