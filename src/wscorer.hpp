#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "forest.hpp"
#include "scorer.hpp"
#include "thread_pool.hpp"
#include "units.hpp"

namespace wqfm {

// Weighted quartets (prototype). The published score (satisfied minus violated quartets, over the
// quartets with two units on each side and no two taxa under one dummy taxon) with every resolved
// quartet weighted as in weighted ASTRAL: by 1 - prod(1 - s_e) over the branches e between its two
// anchor nodes (support weighting), by exp(-L) with L the total length of the four paths from its
// taxa to their anchors (length weighting), or by both. A branch without support counts as fully
// supported, a branch without length as length 0. The real taxa under a dummy taxon keep their Eq. 1
// weights, and the forest is the compressed one of Scorer (built in its forest-only mode).
//
// Every resolved quartet is counted once, at the lowest common ancestor m of its two anchor nodes.
// Either the two pairs of its displayed topology lie in two different child branches of m (a product
// of two pair sums, each discounted by the support along the path from m to the pair's anchor), or
// one pair lies in a child branch and the other pair is split at m, its two taxa in two other branches
// of m (single-taxon sums; one of the two may be the top branch). So a node needs, per child branch,
// the pair sums anchored below it and, for the top branch, only single-taxon sums: everything is
// bottom-up except the discounted single-taxon weight of the top branch. Length weighting enters
// through the discounted weights (a taxon counts exp(-length) along its path to the anchor), which
// propagate up the tree by multiplication.
//
// Gains: a node's contribution is multilinear in the values of each of its branches, so moving a real
// taxon changes it by (gradient with respect to its branch's values) times (the change of those
// values), exactly. The changes are sums over the taxon's ancestors w of its weight discounted to w
// times a value at w (the other units joining it there), with support discounts for the K sums; the
// gradients are summed top-down over the ancestors first (R), so every real taxon's gain is one
// discounted path sum (H), computed for all taxa in one top-down pass. Nodes off the path see the
// taxon in their top branch: a gradient times the discounted weight, summed over subtrees (GA, GB).
// Dummy taxa: every node without the dummy's weight below it sees it only in its top branch, where
// the contribution is linear, so a subtree hanging off the dummy's ancestors changes by the dummy's
// discounted weight there times the subtree's gradient sum. At the ancestors themselves the stored
// values are adjusted for the switch: the pair sums involving the dummy (with the other dummies, W;
// with the units on its own side, SS) are side-independent and computed along its ancestor path, and
// of the terms over pairs of dummies only those with the switched dummy are replaced. Cost per dummy
// and ancestor: O(D L), against O(D L^2) for recomputing the node.
class WScorer : public IScorer {
public:
    struct Options {
        bool support = false, length = false;
        double smin = 0, smax = 1;   // supports are rescaled to (s - smin) / (smax - smin), clipped to [0, 1]
    };
    WScorer(const Forest& f, bool compress, const Options& o);

    void setup(const Subproblem& sp, const DummyPool& pool) override;
    void setup_child(const Subproblem& child, const IScorer& parent, int side) override;
    void set_partition(const std::vector<int8_t>& part_r, const std::vector<int8_t>& part_d) override;
    double score() override;
    double gains(std::vector<double>& real_gain, std::vector<double>& dummy_gain) override;
    void scores(const ScoreBatch& batch, std::vector<double>& out) override;
    void swap_real(int i) override { part_r_[i] = static_cast<int8_t>(1 - part_r_[i]); }
    void swap_dummy(int j) override { part_d_[j] = static_cast<int8_t>(1 - part_d_[j]); }
    void set_pool(ThreadPool* pool) override;
    int part_real(int i) const override { return part_r_[i]; }
    int part_dummy(int j) const override { return part_d_[j]; }
    const SubForest& forest() const override { return fb_.forest(); }
    double scale() const override { return scale_; }
    void set_kappa(double kappa) override;
    double kappa() const override { return 1.0; }
    void set_dummy_gamma(double gamma) override;
    const int8_t* part_r_data() const override { return part_r_.data(); }
    const int8_t* part_d_data() const override { return part_d_.data(); }
    int reals() const override { return R_; }
    int dummies() const override { return d_; }

private:
    // Bottom-up values of one node: the discounted weights of the real taxa (rA, rB) and of all units
    // (a, b) on each side below it; the valid pair sums anchored at or below it, plain (PA, PB, M) and
    // discounted by the support path from the node to the pair's anchor (KA, KB, KM), for A-pairs,
    // B-pairs and mixed pairs; and the mixed pairs grouped by the dummy taxon of their member on that
    // dummy's side (Mk, per entry of the node's dummy list) and by both members' dummies (Mkk, for the
    // dummies p on side A and q on side B of the list: entry (rank of p among the A dummies) * (number of
    // B dummies) + rank of q among the B dummies; at most L*L/4 entries), which the exclusions need.
    // The ranks are those of the bipartition of the pass that wrote the values.
    // The scalars are in NodeVals; Mk, KMk (per entry of the node's dummy list) in mk_, kmk_ at dl_off_,
    // and Mkk, KMkk (the A x B block) in mkk_, kmkk_ at pk_off_, flat so that nodes cost no vectors.
    struct NodeVals {
        double rA = 0, rB = 0, a = 0, b = 0, PA = 0, KA = 0, PB = 0, KB = 0, M = 0, KM = 0;
    };
    struct Scratch {                 // per thread: the branch values of one node, and walk state
        int D = 0, L = 0;
        std::vector<double> rA, rB, a, b, PA, KA, PB, KB, M, KM;   // per branch (D children, then the top)
        std::vector<double> x, Mk, KMk;                            // per dummy and branch
        // per pair of dummies p on side A and q on side B, kPair values summed over the child branches that
        // hold both: their products of branch weights, and the pair sums of (p, q) of those branches (see
        // branch_values), at (rk[p] * Lb + rk[q]) * kPair; the dense pair sums per branch would cost L*L*D
        std::vector<double> pacc;
        std::vector<int> rk, cA, cB;                               // rank of each dummy among its side; a child's A / B entries
        std::vector<int> rc0, rc1, posA, posB;                     // binary_node: ranks in the children, A / B positions
        int La = 0, Lb = 0;                                        // dummies on side A / B in the node's list
        std::vector<double> XC;                                    // per dummy: its weight over the child branches
        std::vector<double> Xall, Sxb, Sxy;                        // node_contrib's sums over all branches
        std::vector<double> An, Bn, Anp, E, SMk, SKMk;
        std::vector<double> gA, gB, gPA, gKA, gPB, gKB, gM, gKM, gMk, gKMk;   // gradients per branch
        std::vector<int> pos;
        std::vector<int8_t> pd;
        // dummy gains: overriding values for the nodes whose values change (ovr: node -> index in tmp,
        // or -1), and an overriding top branch (useTop)
        std::vector<int> ovr, path;
        std::vector<NodeVals> tmp;
        bool useTop = false;
        double aT = 0, bT = 0;
        double pairPart = 0;                 // node_contrib: the terms over pairs of dummies
        // dummy gains: per ancestor of the switched dummy (index in ovr), its pair sums with the other
        // dummies of the node's list (W, KW at woff) and with the other units of its own side (SS, KSS)
        std::vector<double> W, KW;
        struct PathVals { int woff; double SS, KSS; };
        std::vector<PathVals> pv;
        std::vector<double> Wc, KWc;         // at a node: W per branch (row-major [row][branch])
        std::vector<double> Mo, KMo;         // at a node: the switched dummy's stored mixed-pair sums, per row and branch
        std::vector<double> SSc, KSSc, OPc, KOPc;   // per branch
        std::vector<double> bx;              // dummy_change at a binary node: per-row arrays (kRows blocks of L)
        std::vector<int8_t> sd;              // per row: side after the switch
    };

    void build();
    void branch_values(int m, const int8_t* pd, Scratch& S) const;
    void node_values(int m, const int8_t* pd, Scratch& S, NodeVals& out);
    double node_contrib(int m, const int8_t* pd, Scratch& S, bool grad, bool pairs = true) const;
    double binary_node(int v, const int8_t* pd, bool keep, Scratch& S);
    double pair_term(int p, int q, const double* M, const double* K, Scratch& S) const;
    double tree_pass(int g, const int8_t* pr, const int8_t* pd, bool keep, Scratch& S);
    void gain_pass(int g, const int8_t* pr, const int8_t* pd, double* real_gain, Scratch& S);
    double dummy_change(int g, int tq, int j, Scratch& S);
    double dummy_change_binary(int p, int j, int tk, int from, double xt, Scratch& S, Scratch::PathVals& me);
    double dummy_change_poly(int p, int j, int tk, int from, double xt, Scratch& S, Scratch::PathVals& me);
    void run_chunks(const std::function<void(int, int)>& fn);
    void ensure_scratch();

    Scorer fb_;                       // builds the forests
    Options opt_;
    ThreadPool* pool_ = nullptr;
    const SubForest* C_ = nullptr;
    int R_ = 0, d_ = 0;
    std::vector<int8_t> part_r_, part_d_;
    double scale_ = 0;

    // per node, fixed within a subproblem: the length factor exp(-length) and support factor 1 - s of
    // its branch; the dummies with weight below it (ascending), that weight discounted to the node, and
    // their positions among the tree's dummies; per tree its dummies; per internal node the discounted
    // weight of each of the tree's dummies in its top branch
    std::vector<double> ef_, df_;
    std::vector<int> cpos_;                    // per node: its index among its parent's children
    std::vector<int> dl_off_, dl_k_, dl_t_;
    std::vector<double> dl_x_;
    std::vector<int> pres_off_, pres_;
    std::vector<int> cl_off_, cl_;            // per (tree, dummy) entry of pres_: its cluster leaves
    std::vector<int> xtop_off_;
    std::vector<double> xtop_;
    // per node, per pass
    std::vector<NodeVals> nv_;
    std::vector<double> mk_, kmk_, mkk_, kmkk_;          // see NodeVals
    std::vector<size_t> pk_off_;                        // per node: offset of its A x B block (L*L/4 entries)
    std::vector<double> rAtop_, rBtop_, aTop_, bTop_;   // discounted weights in the top branch
    std::vector<double> c_, pp_;                        // contribution, and its terms over pairs of dummies
    std::vector<double> ga_, gb_;                       // the contribution's gradients w.r.t. the top sums
    std::vector<double> GA_, GB_;                       // their sums over subtrees, discounted
    // gains: per node (as a child branch of its parent) the gradients of the parent's contribution
    // w.r.t. its branch values (8: a b PA KA PB KB M KM) and, per entry of the parent's dummy list and
    // child, w.r.t. Mk and KMk; the gradient sums over the ancestors (R*), and the path sums (H)
    std::vector<double> ge_, gmk_, gkmk_;
    std::vector<int> gmk_off_;
    std::vector<double> rPA_, rKA_, rPB_, rKB_, rM_, rKM_, rMk_, rKMk_, H_;
    // chunks of trees
    std::vector<int> chunk_;
    std::vector<double> chunk_score_, chunk_gain_, batch_score_;
    std::vector<Scratch> scratch_;
};

}  // namespace wqfm
