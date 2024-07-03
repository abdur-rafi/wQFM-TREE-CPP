#pragma once

enum ScoreNormalizationType{
    NO_NORMALIZATION,
    FLAT_NORMALIZATION,
    NESTED_NORMALIZATION
};

enum ConsensusWeightType{
    FLAT,
    NESTED
};

const int N_THREADS = 1;

const ConsensusWeightType CONSENSUS_WEIGHT_TYPE = ConsensusWeightType::NESTED;

double scoreEqn(double total, double sat){
    return 2 * sat - total;
}

const bool USE_SCORING_IN_CONSENSUS = true;

const int MAX_ITERATION = 5;

const bool DEBUG = false;