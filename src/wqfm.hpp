#pragma once

#include <memory>
#include <string>
#include <vector>

#include "consensus.hpp"
#include "forest.hpp"
#include "scorer.hpp"
#include "trace.hpp"
#include "units.hpp"

namespace wqfm {

// Unrooted species tree (or subtree with dummy-taxon leaves) as an adjacency list.
struct SpeciesTree {
    std::vector<std::vector<int>> adj;
    std::vector<int> taxon;    // real taxon at a leaf, else -1
    std::vector<int> dummy;    // dummy taxon at a leaf, else -1

    int add(int t, int d);
    std::string newick(const Taxa& taxa) const;
};

struct SearchStats {
    long subproblems = 0, passes = 0, moves = 0;
    double t_setup = 0, t_initial = 0, t_fm = 0;   // seconds: scorer setup, initial bipartitions, FM passes
};

// The wQFM divide-and-conquer (paper Section 2, Supplementary Fig. 1), replicating QFM2.java:
// consensus-based initial bipartition, FM passes (at most MAX_ITERATION + 1 = 6) with the Java
// tie-breaking, divide with one new dummy taxon per side, and joining the two subtrees.
class Search {
public:
    Search(const Forest& f, const ConsensusPartitioner& cons, Trace* trace, const Replay* replay, bool compress = true)
        : F_(f), cons_(cons), trace_(trace), replay_(replay), compress_(compress) {}

    SpeciesTree run();
    const SearchStats& stats() const { return stats_; }

private:
    // a subproblem's scorer is set up once, used for its initial bipartition and then for its FM passes
    // (child subproblems derive their forest from the parent's: pass the parent scorer and the side)
    std::unique_ptr<Scorer> prepare(Subproblem& sp, const Scorer* parent = nullptr, int side = 0);
    SpeciesTree recurse(Subproblem& sp, std::unique_ptr<Scorer> sc);
    bool fm_pass(Subproblem& sp, Scorer& sc, const Replay::Pass* rp);

    const Forest& F_;
    const ConsensusPartitioner& cons_;
    Trace* trace_;
    const Replay* replay_;
    bool compress_;
    std::vector<std::unique_ptr<Scorer>> spare_;   // reused to keep their buffers
    DummyPool pool_;
    SearchStats stats_;
};

}  // namespace wqfm
