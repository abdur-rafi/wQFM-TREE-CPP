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

// The gene trees as seen by one subproblem. Leaves are real taxa or "clusters": a dummy taxon with the
// summed weight of its taxa there. Children precede parents within a tree and the root is last.
//
// With compression, every maximal part of a gene tree whose taxa all lie under one dummy taxon becomes
// one cluster (judged on the unrooted tree, from below and above every node), clusters of the same
// dummy hanging off one node are merged, and gene trees with fewer than 4 units are dropped. All of
// this is exact: the removed nodes only anchor quartets with two taxa under one dummy, which the score
// excludes, and the dropped trees contribute nothing to any score or gain (tools CompressCheck and
// brute_check verify it).
struct SubForest {
    int k = 0;
    std::vector<int> tree_off, parent, child_off, children;
    std::vector<int> leaf_real;      // real unit index at a leaf, else -1
    std::vector<int> leaf_dummy;     // dummy index at a cluster leaf, else -1
    std::vector<double> leaf_w;      // leaf weight (1 for a real taxon)
    std::vector<int> real_leaf;      // k*R: node of real unit i in tree g, or -1

    int nodes() const { return static_cast<int>(parent.size()); }
};

// Score of a bipartition and the gain of moving each unit to the other side, computed directly from
// the gene trees (paper Section 2.5).
class Scorer {
public:
    explicit Scorer(const Forest& f, bool compress = true) : F_(f), compress_(compress) {}

    void setup(const Subproblem& sp, const DummyPool& pool);   // from the gene trees
    // from the parent's forest: `child` is side `side` of the parent's final bipartition plus the new dummy
    void setup_child(const Subproblem& child, const Scorer& parent, int side);
    void set_partition(const std::vector<int8_t>& part_r, const std::vector<int8_t>& part_d);

    double score();
    // real_gain[i] (dummy_gain[j]) = score change if real taxon i (dummy j) moves to the other side
    double gains(std::vector<double>& real_gain, std::vector<double>& dummy_gain);

    void swap_real(int i);
    void swap_dummy(int j);

    int part_real(int i) const { return part_r_[i]; }
    int part_dummy(int j) const { return part_d_[j]; }
    const SubForest& forest() const { return C_; }

private:
    void begin_forest();
    void finish_setup();
    template <class Src> void copy_from(const Src& src, int g);
    template <class Src> void compress_from(const Src& src, int g);
    double tree_pass(int g, bool want, double* real_gain, double* dummy_gain);
    double generic_value(int D, int moved, PolyEval& pe, const std::vector<int>& pres);

    const Forest& F_;
    bool compress_;
    int d_ = 0, R_ = 0;
    std::vector<int> unit_r_, unit_d_;        // per taxon: real index / dummy index, or -1
    std::vector<double> w_;                   // per taxon weight (Eq. 1)
    std::vector<int8_t> part_r_, part_d_;
    SubForest C_;

    // per node of C_: below-subtree totals
    std::vector<double> a_, b_;               // weight in partition 0 / 1
    std::vector<int> ra_, rb_;                // real taxa in partition 0 / 1
    std::vector<double> x_;                   // nodes*d weight of each dummy (fixed within a subproblem)
    // per tree of C_
    std::vector<double> ta_, tb_, tx_;
    std::vector<int> tra_, trb_;
    std::vector<std::vector<int>> present_;   // dummies with weight in the tree

    // scratch
    std::vector<double> push_, dsum_, xt_;
    std::vector<double> ga_, gb_;             // generic node branch values
    std::vector<int> gra_, grb_;
    std::vector<const double*> gx_;
    std::vector<double> gd_;
    PolyEval pe_, pe2_;
    std::vector<double> agg_backup_;
    // compression scratch
    std::vector<int> lab_down_, lab_up_, cid_, stamp_, slot_, suffix_, touched_, ni_;
    std::vector<std::pair<int, double>> clusters_;
    std::vector<double> wsum_;
    std::vector<char> kept_;
    std::vector<std::vector<int>> adj_;
    std::vector<int> node_real_, node_dummy_;
    std::vector<double> node_w_;
};

}  // namespace wqfm
