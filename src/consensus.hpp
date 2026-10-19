#pragma once

#include <vector>

#include "forest.hpp"
#include "scorer.hpp"
#include "units.hpp"

namespace wqfm {

class Trace;
struct Replay;

// Initial bipartition of a subproblem from a greedy consensus tree (paper Section 2.4, Supplementary
// Algorithm 1), replicating the Java ConsensusTreePartition: every consensus edge that splits the units
// into two sides of at least 2 (or 1 with allow_singleton) is scored, and the first best one wins.
// The candidate order decides ties, so the tree keeps the Java parser's reversed child order and the
// Java root re-balancing.
class ConsensusPartitioner {
public:
    ConsensusPartitioner(const ParsedTree& tree, const Taxa& taxa);   // throws if a taxon is missing

    // candidates are compared by (2F - kappa T) / T^power (Search options); default: the scorer's score
    void set_objective(double kappa, double power) { kappa_ = kappa; power_ = power; }
    // score the candidates with this kappa instead of the scorer's (linear objectives; < 0: off)
    void set_init_kappa(double kappa) { init_kappa_ = kappa; }

    // sets sp.part_r and sp.part_d; `scorer` must be set up for sp
    void make(Subproblem& sp, bool allow_singleton, IScorer& scorer, const DummyPool& pool,
              Trace* trace, const Replay* replay) const;

private:
    void balance_root();

    std::vector<std::vector<int>> childs_;
    std::vector<int> parent_, taxon_, topo_;
    int root_ = -1;
    int n_ = 0;
    double kappa_ = 1, power_ = 0, init_kappa_ = -1;
};

}  // namespace wqfm
