#pragma once

enum ScoreNormalizationType{
    NO_NORMALIZATION,
    FLAT_NORMALIZATION,
    NESTED_NORMALIZATION
};

const int N_THREADS = 1;


double scoreEqn(double total, double sat){
    return 2 * sat - total;
}