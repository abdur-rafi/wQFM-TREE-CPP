#pragma once

#include <cstdint>
#include <vector>

#include "forest.hpp"
#include "units.hpp"

namespace wqfm {

// Power-sum evaluation of one node's contribution f = w(S(g,u)) + w(U(g,u))/2 (paper Sections 2.5.1
// and 2.5.3; the Java "sat" per node). Each branch i carries component values: alpha_0 = weight in
// partition A, alpha_m = weight of the m-th dummy in A (sign -1), and likewise beta for B, plus rho =
// real-taxon count in B. Then PA_ij = sum_m s_m alpha_mi alpha_mj, PB_kl likewise, and
// PB_k = (sum_n s_n beta_nk^2 - rho_k)/2, so f is a function of additive per-branch aggregates.
struct PolyEval {
    int D = 0, mA = 0, mB = 0;
    std::vector<double> al, be, rho;   // al[m*D+i], be[n*D+i], rho[i]
    std::vector<double> agg;           // TA QA R R2 (mA each), TB QB (mB each), M1..M4 (mA*mB each), P

    void resize(int d, int ma, int mb);
    void build();                      // aggregates from al/be/rho
    void add(int i, double s);         // add (s=+1) or remove (s=-1) branch i
    double eval() const;
};

// Score of a bipartition and the gain of moving each unit to the other side, computed directly from
// the gene trees. Phase 1: exact and straightforward (dense per-dummy weights, no compression).
class Scorer {
public:
    explicit Scorer(const Forest& f) : F_(f) {}

    void setup(const Subproblem& sp, const DummyPool& pool);   // units, weights, per-dummy weights; O(N d)
    void set_partition(const std::vector<int8_t>& part_r, const std::vector<int8_t>& part_d);   // O(N)

    double score();
    // real_gain[i] (dummy_gain[j]) = score change if real taxon i (dummy j) moves to the other side
    double gains(std::vector<double>& real_gain, std::vector<double>& dummy_gain);

    void swap_real(int i);
    void swap_dummy(int j);

    int part_real(int i) const { return part_r_[i]; }
    int part_dummy(int j) const { return part_d_[j]; }

private:
    double tree_pass(int g, bool want, double* real_gain, double* dummy_gain);
    double generic_value(int D, int moved, PolyEval& pe, const std::vector<int>& pres);

    const Forest& F_;
    int d_ = 0;
    std::vector<int> reals_;                  // taxon id of each real unit
    std::vector<int> unit_r_, unit_d_;        // per taxon: real index / dummy index, or -1
    std::vector<double> w_;                   // per taxon weight (Eq. 1)
    std::vector<int8_t> part_r_, part_d_;

    // per node: below-subtree totals
    std::vector<double> a_, b_;               // weight in partition 0 / 1
    std::vector<int> ra_, rb_;                // real taxa in partition 0 / 1
    std::vector<double> x_;                   // N*d weight of each dummy (fixed within a subproblem)
    // per tree
    std::vector<double> ta_, tb_, tx_;
    std::vector<int> tra_, trb_;
    std::vector<std::vector<int>> present_;   // dummies with weight in the tree

    // scratch
    std::vector<double> push_, dsum_, xt_;
    std::vector<double> ga_, gb_;             // generic node branch values
    std::vector<int> gra_, grb_;
    std::vector<const double*> gx_;
    std::vector<double> gd_;                  // generic node deltas (2 per branch)
    PolyEval pe_, pe2_;
    std::vector<double> agg_backup_;
};

}  // namespace wqfm
