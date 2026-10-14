#pragma once

#include <cstddef>
#include <cstdint>
#include <new>
#include <vector>

#include "forest.hpp"
#include "thread_pool.hpp"
#include "units.hpp"

namespace wqfm {

// Allocator for per-thread scratch: cache-line aligned and padded, so buffers written by different
// threads never share a cache line (small adjacent heap blocks would, and false sharing then
// dominates the per-node work).
template <class T>
struct LineAlloc {
    using value_type = T;
    static constexpr std::size_t kLine = 64;
    LineAlloc() = default;
    template <class U> LineAlloc(const LineAlloc<U>&) {}
    T* allocate(std::size_t n) {
        return static_cast<T*>(::operator new((n * sizeof(T) + kLine - 1) / kLine * kLine, std::align_val_t(kLine)));
    }
    void deallocate(T* p, std::size_t) { ::operator delete(p, std::align_val_t(kLine)); }
    template <class U> bool operator==(const LineAlloc<U>&) const { return true; }
};
template <class T> using LineVec = std::vector<T, LineAlloc<T>>;

// Contribution f = w(S(g,u)) + w(U(g,u))/2 of one polytomy node u in a tree with dummy taxa (paper
// Sections 2.5.1 and 2.5.3; the Java "sat" per node), and its change when one real taxon or one dummy
// taxon switches sides.
//
// Branch i carries a_i and b_i (weight in partitions A and B), rho_i (real taxa in B) and x_ki (weight
// of present dummy k). With PA_ij = a_i a_j - sum_{k in A} x_ki x_kj, PB_i = (b_i^2 - rho_i - sum_{k in B}
// x_ki^2)/2 and PB_ij = b_i b_j - sum_{k in B} x_ki x_kj,
//   S = sum_{i<j} PA_ij sum_{l != i,j} PB_l,   U = sum_{i<j} PA_ij sum_{l<m; l,m != i,j} PB_lm,
// which expand into power sums over the branches: of the "components" (all of A, all of B, each dummy)
// and of pairs of components on opposite sides. A dummy's weight lies in a few branches, so the sums
// are kept per dummy and per pair of dummies, together with their totals over each side. Moving a real
// taxon then costs O(1 + dummies in its branch) and moving a dummy O(d + its branches), instead of
// rebuilding everything (O(D d^2)) for each move.
struct PolyNode {
    // input: begin(), then branch() for each branch followed by entry() for its nonzero dummy weights
    void begin(int dummies);
    void branch(double a, double b, double rho, double sq_extra = 0);   // sq_extra: B-side squares of
                                                                       // dummies not given as entries
    void entry(int k, double x);
    LineVec<int> inA;                          // per dummy: 1 if on side A (set before build)

    double build();                            // aggregates; returns f
    double real_delta(int i, int dir);         // change of f when a real taxon in branch i moves A->B (0) or B->A (1)
    double dummy_delta(int k);                 // change of f when dummy k switches sides
    double pairs_without(int i) const;         // sum of PA over the pairs of branches other than i

private:
    struct Sums {
        double TA0, QA0, R0, R20, TB0, QB0, P, N1, N2, N3, N4;   // components "all of A" and "all of B"
        double EAd, EBd, srd, CG1, CG0, CH1, CH0, Xp;            // dummy terms summed over each side
    };
    static constexpr int kAgg = 10;            // per dummy: R R2 G1..G4 H1..H4
    static double eval(const Sums& s);
    static void terms(Sums& s, bool in_a, double t, double q, const double* g, double sign);
    void branch_change(Sums& s, int i, double a2, double b2, double rho2, double sqB2) const;

    int p_ = 0;
    LineVec<double> a_, b_, rho_, sq0_, sqB_, pb_;   // per branch
    LineVec<int> eoff_, ek_;                   // entries of branch i: [eoff_[i], eoff_[i+1])
    LineVec<double> ex_;
    LineVec<int> koff_, kbr_, kfill_;          // per dummy: its branches and weights
    LineVec<double> kx_;
    LineVec<double> T_, Q_, agg_, pair_;       // pair_[(k*p+l)*4]: sums of x_k x_l, x_k^2 x_l, x_k x_l^2, x_k^2 x_l^2
    LineVec<double> pairA_, pairB_;            // per dummy: its opposite-side pair terms as a member of A / of B
    Sums s_{};
    double f0_ = 0;
    LineVec<double> dacc_;                     // dummy_delta: per dummy, change of its sums (kAgg each)
    LineVec<int> dstamp_, dtouched_;
    int tag_ = 0;
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
//
// The trees are processed in chunks fixed by the forest alone (never by the thread count) and per-chunk
// results are summed in chunk order, so results are identical with or without a thread pool.
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

    // move a unit to the other side; the totals follow in the next score() or gains(), chunk by chunk
    void swap_real(int i);
    void swap_dummy(int j);

    // run the per-tree work on this pool (nullptr: on the calling thread only)
    void set_pool(ThreadPool* pool);

    int part_real(int i) const { return part_r_[i]; }
    int part_dummy(int j) const { return part_d_[j]; }
    const SubForest& forest() const { return C_; }
    // sum over trees of W^4, W = the tree's total unit weight: the scale of the quartet weights a score
    // sums, hence of its rounding error
    double scale() const { return scale_; }

private:
    void begin_forest();
    void finish_setup();
    template <class Src> void copy_from(const Src& src, int g);
    template <class Src> void compress_from(const Src& src, int g);
    struct Scratch {                          // per thread
        LineVec<double> dsum, gd, pdelta;          // dsum, pos: per dummy
        LineVec<int> pos;
        PolyNode pn;
    };

    double tree_pass(int g, bool want, double* real_gain, double* dummy_gain, Scratch& S);
    double tree_pass_plain(int g, bool want, double* real_gain);
    double tree_pass_dummies(int g, bool want, double* real_gain, double* dummy_gain, Scratch& S);
    double poly_node(int v, const double* X, double sqXB, double at, double bt, int rat, int rbt, bool want,
                     double& G0, double& G1, double coef[3], Scratch& S);
    void run_chunks(const std::function<void(int, int)>& fn);
    void ensure_scratch();
    void set_partition_chunk(int c);
    void apply_moves(int c);

    const Forest& F_;
    bool compress_;
    ThreadPool* pool_ = nullptr;
    int d_ = 0, R_ = 0;
    std::vector<int> unit_r_, unit_d_;        // per taxon: real index / dummy index, or -1
    std::vector<double> w_;                   // per taxon weight (Eq. 1)
    std::vector<int8_t> part_r_, part_d_;
    // moves not yet applied to the totals below, in order (part_r_ / part_d_ already include them)
    struct Move { int unit; bool dummy; int8_t from; };
    std::vector<Move> moves_;
    SubForest C_;
    double scale_ = 0;

    // per node of C_: below-subtree totals
    std::vector<double> a_, b_;               // weight in partition 0 / 1
    std::vector<int> ra_, rb_;                // real taxa in partition 0 / 1
    // per node: the dummies with weight below it (ascending) and that weight, and at a binary node also
    // their weights below its two children; per dummy: its nodes (ascending) and weights. Fixed within a
    // subproblem.
    std::vector<int> xoff_, xj_, doff_, dnode_;
    std::vector<double> xw_, xc0_, xc1_, dw_;
    // per tree of C_
    std::vector<double> ta_, tb_, tx_;
    std::vector<int> tra_, trb_;
    std::vector<std::vector<int>> present_;   // dummies with weight in the tree

    // per tree (offsets inner_off_ / rleaf_off_): its internal nodes in post-order and its real-taxon
    // leaves, so the kernels loop without testing for leaves
    std::vector<int> inner_, inner_off_, rleaf_, rleaf_off_;
    // per node: values pushed down to compute real-taxon gains
    std::vector<double> push_;
    // chunks of trees (fixed by the forest) and their partial results
    std::vector<int> chunk_;
    std::vector<double> chunk_score_, chunk_gain_;
    std::vector<Scratch> scratch_;
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
