#pragma once

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "consensus.hpp"
#include "forest.hpp"
#include "scorer.hpp"
#include "thread_pool.hpp"
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

struct SearchOptions {
    bool compress = true;      // compressed subproblem forests (exact; off only for testing)
    int threads = 1;
    // FM keeps a longer prefix of moves only if its cumulative gain beats the best so far by more than
    // 1e-5 (Java's EPS). With relative_eps the margin is 1e-5 + 1e-12 * (sum of |gain| of the moves so
    // far): at large score scales, rounding in the running sum can exceed 1e-5 (in Java a pass that
    // moved every unit, i.e. the same bipartition, summed to +2.9e-4 and counted as an improvement).
    bool relative_eps = false;
};

struct SearchStats {
    long subproblems = 0, passes = 0, moves = 0;
    double t_setup = 0, t_initial = 0, t_fm = 0;   // seconds summed over threads
    int threads = 1;
};

// The wQFM divide-and-conquer (paper Section 2, Supplementary Fig. 1), replicating QFM2.java:
// consensus-based initial bipartition, FM passes (at most MAX_ITERATION + 1 = 6) with the Java
// tie-breaking, divide with one new dummy taxon per side, and joining the two subtrees.
//
// With several threads, each subproblem spreads its per-tree work over the pool, and the two children
// of a divide step run concurrently (one inline, one as a task). Every subproblem is computed the same
// way whatever thread runs it, so the tree does not depend on the thread count. Tracing and replay
// process subproblems one at a time (their per-tree work is still parallel).
class Search {
public:
    Search(const Forest& f, const ConsensusPartitioner& cons, Trace* trace, const Replay* replay,
           const SearchOptions& opt = {});

    SpeciesTree run();
    SearchStats stats() const;

private:
    // a subproblem's scorer is set up once, used for its initial bipartition and then for its FM passes
    // (child subproblems derive their forest from the parent's: pass the parent scorer and the side)
    std::unique_ptr<Scorer> prepare(Subproblem& sp, const Scorer* parent, int side);
    SpeciesTree recurse(Subproblem& sp, std::unique_ptr<Scorer> sc);
    bool fm_pass(Subproblem& sp, Scorer& sc, const Replay::Pass* rp);
    std::unique_ptr<Scorer> acquire();
    void release(std::unique_ptr<Scorer> sc);
    void add_time(double& total, double secs);

    const Forest& F_;
    const ConsensusPartitioner& cons_;
    Trace* trace_;
    const Replay* replay_;
    bool compress_, relative_eps_;
    std::unique_ptr<ThreadPool> threads_;
    bool tasks_ = false;
    DummyPool pool_;
    mutable std::mutex mu_;                       // spare_ and the times
    std::vector<std::unique_ptr<Scorer>> spare_;  // reused to keep their buffers
    std::atomic<long> subproblems_{0}, passes_{0}, moves_{0};
    double t_setup_ = 0, t_initial_ = 0, t_fm_ = 0;
};

}  // namespace wqfm
