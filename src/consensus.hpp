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

    // sets sp.part_r and sp.part_d
    void make(Subproblem& sp, bool allow_singleton, Scorer& scorer, const DummyPool& pool,
              Trace* trace, const Replay* replay) const;

private:
    void balance_root();

    std::vector<std::vector<int>> childs_;
    std::vector<int> parent_, taxon_, topo_;
    int root_ = -1;
    int n_ = 0;
};

}  // namespace wqfm
