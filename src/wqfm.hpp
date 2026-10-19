#pragma once

#include <atomic>
#include <chrono>
#include <cstdio>
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
    double progress_secs = 0;  // print the search progress to stderr this often (0: never)
    // FM keeps a longer prefix of moves only if its cumulative gain beats the best so far by more than
    // 1e-5 (Java's EPS). With relative_eps the margin is 1e-5 + 1e-12 * (sum of |gain| of the moves so
    // far): at large score scales, rounding in the running sum can exceed 1e-5 (in Java a pass that
    // moved every unit, i.e. the same bipartition, summed to +2.9e-4 and counted as an improvement).
    bool relative_eps = false;
    // if set, one line per subproblem: real taxa, dummy taxa, score of the initial bipartition, score
    // after the FM passes, FM passes, and 1 if FM changed the bipartition (lines in completion order)
    std::FILE* subproblem_log = nullptr;
    // Experimental objectives (defaults: the published wQFM-TREE): FM and the initial bipartition
    // maximise (2F - kappa T) / T^power per subproblem, with F and T as in Scorer::set_kappa and
    // dummy-taxon weights by dummy_gamma (Scorer::set_dummy_gamma).
    double kappa = 1, power = 0, dummy_gamma = 0;
    double init_kappa = -1;    // kappa for scoring the initial-bipartition candidates (< 0: kappa)
    // weighted quartets (WScorer): 0 off, 1 support, 2 length, 3 both; supports rescaled from [smin, smax]
    int weighting = 0;
    double support_min = 0, support_max = 1;
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
    std::unique_ptr<IScorer> prepare(Subproblem& sp, const IScorer* parent, int side);
    SpeciesTree recurse(Subproblem& sp, std::unique_ptr<IScorer> sc);
    bool fm_pass(Subproblem& sp, IScorer& sc, const Replay::Pass* rp);
    std::unique_ptr<IScorer> acquire();
    void release(std::unique_ptr<IScorer> sc);
    void add_time(double& total, double secs);
    void report_progress();

    const Forest& F_;
    const ConsensusPartitioner& cons_;
    Trace* trace_;
    const Replay* replay_;
    bool compress_, relative_eps_;
    double progress_secs_ = 0;
    std::chrono::steady_clock::time_point t_start_;
    std::atomic<long long> next_report_ns_{0};   // since t_start_
    std::unique_ptr<ThreadPool> threads_;
    bool tasks_ = false;
    DummyPool pool_;
    mutable std::mutex mu_;                       // the times
    std::FILE* sublog_ = nullptr;
    std::mutex log_mu_;
    double kappa_ = 1, power_ = 0, gamma_ = 0;
    int weighting_ = 0;
    double smin_ = 0, smax_ = 1;
    double objective_gains(IScorer& sc, std::vector<double>& gr, std::vector<double>& gd);
    double objective_value(IScorer& sc);
    std::atomic<long> subproblems_{0}, passes_{0}, moves_{0};
    double t_setup_ = 0, t_initial_ = 0, t_fm_ = 0;
};

}  // namespace wqfm
