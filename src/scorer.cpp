#include "scorer.hpp"

#include <algorithm>
#include <stdexcept>

namespace wqfm {

// ---------------------------------------------------------------- PolyEval
//
// With s_m = +1 for m = 0 and -1 for dummy components:
//   EA = sum_{i<j} PA_ij = sum_m s_m (TA_m^2 - QA_m)/2,  EB likewise
//   S  = sum_{i<j} PA_ij (P - PB_i - PB_j)               = P EA - sum_m s_m (TA_m R_m - R2_m)
//   U  = sum_{i<j} PA_ij sum_{k<l, k,l not in {i,j}} PB_kl = EA EB - sum_i SB_i g_i + sum_{i<j} PA_ij PB_ij
// where g_i = sum_{j!=i} PA_ij and SB_i = sum_{l!=i} PB_il, which expand into the M1..M4 power sums.

namespace {
inline double sgn(int m) { return m == 0 ? 1.0 : -1.0; }
}

void PolyEval::resize(int d, int ma, int mb) {
    D = d; mA = ma; mB = mb;
    al.assign(static_cast<size_t>(ma) * d, 0.0);
    be.assign(static_cast<size_t>(mb) * d, 0.0);
    rho.assign(d, 0.0);
    agg.assign(4 * ma + 2 * mb + 4 * ma * mb + 1, 0.0);
}

void PolyEval::build() {
    std::fill(agg.begin(), agg.end(), 0.0);
    for (int i = 0; i < D; ++i) add(i, 1.0);
}

void PolyEval::add(int i, double s) {
    double* TA = agg.data();
    double* QA = TA + mA;
    double* R = QA + mA;
    double* R2 = R + mA;
    double* TB = R2 + mA;
    double* QB = TB + mB;
    double* M1 = QB + mB;
    double* M2 = M1 + mA * mB;
    double* M3 = M2 + mA * mB;
    double* M4 = M3 + mA * mB;
    double* P = M4 + mA * mB;

    double pb = -rho[i];
    for (int n = 0; n < mB; ++n) {
        const double b = be[static_cast<size_t>(n) * D + i];
        pb += sgn(n) * b * b;
    }
    pb *= 0.5;
    for (int m = 0; m < mA; ++m) {
        const double a = al[static_cast<size_t>(m) * D + i];
        if (a == 0) continue;
        TA[m] += s * a;
        QA[m] += s * a * a;
        R[m] += s * pb * a;
        R2[m] += s * pb * a * a;
        for (int n = 0; n < mB; ++n) {
            const double b = be[static_cast<size_t>(n) * D + i];
            if (b == 0) continue;
            const int k = m * mB + n;
            M1[k] += s * a * b;
            M2[k] += s * a * a * b;
            M3[k] += s * a * b * b;
            M4[k] += s * a * a * b * b;
        }
    }
    for (int n = 0; n < mB; ++n) {
        const double b = be[static_cast<size_t>(n) * D + i];
        TB[n] += s * b;
        QB[n] += s * b * b;
    }
    *P += s * pb;
}

double PolyEval::eval() const {
    const double* TA = agg.data();
    const double* QA = TA + mA;
    const double* R = QA + mA;
    const double* R2 = R + mA;
    const double* TB = R2 + mA;
    const double* QB = TB + mB;
    const double* M1 = QB + mB;
    const double* M2 = M1 + mA * mB;
    const double* M3 = M2 + mA * mB;
    const double* M4 = M3 + mA * mB;
    const double P = *(M4 + mA * mB);

    double EA = 0, EB = 0, sr = 0, X = 0;
    for (int m = 0; m < mA; ++m) {
        EA += sgn(m) * (TA[m] * TA[m] - QA[m]) / 2;
        sr += sgn(m) * (TA[m] * R[m] - R2[m]);
    }
    for (int n = 0; n < mB; ++n) EB += sgn(n) * (TB[n] * TB[n] - QB[n]) / 2;
    for (int m = 0; m < mA; ++m)
        for (int n = 0; n < mB; ++n) {
            const int k = m * mB + n;
            X += sgn(m) * sgn(n) *
                 (TB[n] * TA[m] * M1[k] - TB[n] * M2[k] - TA[m] * M3[k] + M4[k] - (M1[k] * M1[k] - M4[k]) / 2);
        }
    const double S = P * EA - sr;
    const double U = EA * EB - X;
    return S + U / 2;
}

// ---------------------------------------------------------------- Scorer

namespace {

// Node contribution with no dummy taxa in the tree (all weight real, so rho = b): fixed-size power
// sums, the PolyEval formulas with one component on each side.
struct Agg0 { double A, B, Saa, Sbb, Sp, Sap, Saap, Sab, Sabb, Saab, Saabb; };

inline void agg0_add(Agg0& g, double a, double b, double s) {
    const double p = b * (b - 1) / 2;
    g.A += s * a; g.B += s * b; g.Saa += s * a * a; g.Sbb += s * b * b; g.Sp += s * p;
    g.Sap += s * a * p; g.Saap += s * a * a * p; g.Sab += s * a * b; g.Sabb += s * a * b * b;
    g.Saab += s * a * a * b; g.Saabb += s * a * a * b * b;
}

inline double agg0_f(const Agg0& g) {
    const double E2a = (g.A * g.A - g.Saa) / 2, E2b = (g.B * g.B - g.Sbb) / 2;
    const double S = g.Sp * E2a - g.A * g.Sap + g.Saap;
    const double U = E2a * E2b - (g.A * (g.B * g.Sab - g.Sabb) - (g.B * g.Saab - g.Saabb)) + (g.Sab * g.Sab - g.Saabb) / 2;
    return S + U / 2;
}

constexpr int kMixed = -1, kEmpty = -2;
inline int combine(int a, int b) {
    if (a == kEmpty) return b;
    if (b == kEmpty) return a;
    return (a == b && a >= 0) ? a : kMixed;
}
}  // namespace

// Sources of trees for building a subproblem forest: the gene trees themselves, or the parent
// subproblem's forest relabelled with the child's units. Building from the parent is exact: the
// child's units are unions of the parent's (one side keeps its units, the other side becomes the new
// dummy taxon), so whatever was pure under one parent unit is pure under one child unit, and a tree
// with fewer than 4 parent units has fewer than 4 child units.
struct LeafUnit { int real, dummy; double w; };

struct GeneTreeSource {
    const Forest& F;
    const std::vector<int>& unit_r;
    const std::vector<int>& unit_d;
    const std::vector<double>& w;
    int trees() const { return F.k; }
    int begin(int g) const { return F.tree_off[g]; }
    int end(int g) const { return F.tree_off[g + 1]; }
    int cbeg(int v) const { return F.child_off[v]; }
    int cend(int v) const { return F.child_off[v + 1]; }
    int child(int e) const { return F.children[e]; }
    int parent(int v) const { return F.parent[v]; }
    bool leaf(int v) const { return F.taxon[v] >= 0; }
    LeafUnit unit(int v) const {
        const int t = F.taxon[v];
        return unit_r[t] >= 0 ? LeafUnit{unit_r[t], -1, 1.0} : LeafUnit{-1, unit_d[t], w[t]};
    }
};

struct ParentSource {
    const SubForest& P;
    const std::vector<int>& real_map;    // parent real unit -> child real unit, or -1 (into the new dummy)
    const std::vector<int>& dummy_map;   // parent dummy -> child dummy, or -1 (nested in the new dummy)
    int fresh;                           // index of the new dummy taxon
    double m;                            // its child count: units on the other side (Eq. 1)
    int trees() const { return P.k; }
    int begin(int g) const { return P.tree_off[g]; }
    int end(int g) const { return P.tree_off[g + 1]; }
    int cbeg(int v) const { return P.child_off[v]; }
    int cend(int v) const { return P.child_off[v + 1]; }
    int child(int e) const { return P.children[e]; }
    int parent(int v) const { return P.parent[v]; }
    bool leaf(int v) const { return P.child_off[v] == P.child_off[v + 1]; }
    LeafUnit unit(int v) const {
        if (P.leaf_real[v] >= 0) {
            const int r = real_map[P.leaf_real[v]];
            return r >= 0 ? LeafUnit{r, -1, 1.0} : LeafUnit{-1, fresh, 1.0 / m};
        }
        const int j = dummy_map[P.leaf_dummy[v]];
        return j >= 0 ? LeafUnit{-1, j, P.leaf_w[v]} : LeafUnit{-1, fresh, P.leaf_w[v] / m};
    }
};

void Scorer::setup(const Subproblem& sp, const DummyPool& pool) {
    const int n = F_.n;
    d_ = static_cast<int>(sp.dummies.size());
    R_ = static_cast<int>(sp.reals.size());
    unit_r_.assign(n, -1);
    unit_d_.assign(n, -1);
    std::vector<double> coeff(n, 0.0);
    for (int i = 0; i < R_; ++i) {
        unit_r_[sp.reals[i]] = i;
        coeff[sp.reals[i]] = 1.0;
    }
    for (int j = 0; j < d_; ++j) {
        pool.coeffs(sp.dummies[j], 1.0, coeff);
        for (int t : pool[sp.dummies[j]].flat) unit_d_[t] = j;
    }
    w_.assign(n, 0.0);
    for (int t = 0; t < n; ++t) {
        if (unit_r_[t] < 0 && unit_d_[t] < 0) throw std::logic_error("subproblem does not cover every taxon");
        w_[t] = 1. / coeff[t];
    }
    begin_forest();
    const GeneTreeSource src{F_, unit_r_, unit_d_, w_};
    for (int g = 0; g < src.trees(); ++g) {
        if (compress_) compress_from(src, g);
        else copy_from(src, g);
    }
    finish_setup();
}

void Scorer::setup_child(const Subproblem& child, const Scorer& parent, int side) {
    if (!compress_) throw std::logic_error("setup_child requires compression");
    d_ = static_cast<int>(child.dummies.size());
    R_ = static_cast<int>(child.reals.size());
    std::vector<int> real_map(parent.R_, -1), dummy_map(parent.d_, -1);
    int r = 0, j = 0;
    for (int i = 0; i < parent.R_; ++i) if (parent.part_r_[i] == side) real_map[i] = r++;
    for (int i = 0; i < parent.d_; ++i) if (parent.part_d_[i] == side) dummy_map[i] = j++;
    if (r != R_ || j != d_ - 1) throw std::logic_error("setup_child: child does not match the parent's side");
    const double m = (parent.R_ - r) + (parent.d_ - j);
    begin_forest();
    const ParentSource src{parent.C_, real_map, dummy_map, d_ - 1, m};
    for (int g = 0; g < src.trees(); ++g) compress_from(src, g);
    finish_setup();
}

void Scorer::begin_forest() {
    C_ = SubForest{};
    C_.tree_off.push_back(0);
    C_.child_off.push_back(0);
    stamp_.assign(d_, 0);
    slot_.assign(d_, -1);
}

void Scorer::finish_setup() {
    const int N = C_.nodes(), K = C_.k;
    x_.assign(static_cast<size_t>(N) * d_, 0.0);
    if (d_ > 0) {
        for (int v = 0; v < N; ++v) {
            double* xv = &x_[static_cast<size_t>(v) * d_];
            if (C_.leaf_dummy[v] >= 0) {
                xv[C_.leaf_dummy[v]] = C_.leaf_w[v];
            } else {
                for (int e = C_.child_off[v]; e < C_.child_off[v + 1]; ++e) {
                    const double* xc = &x_[static_cast<size_t>(C_.children[e]) * d_];
                    for (int jj = 0; jj < d_; ++jj) xv[jj] += xc[jj];
                }
            }
        }
    }
    tx_.assign(static_cast<size_t>(K) * d_, 0.0);
    present_.assign(K, {});
    for (int g = 0; g < K; ++g) {
        const int root = C_.tree_off[g + 1] - 1;
        for (int jj = 0; jj < d_; ++jj) {
            tx_[static_cast<size_t>(g) * d_ + jj] = x_[static_cast<size_t>(root) * d_ + jj];
            if (tx_[static_cast<size_t>(g) * d_ + jj] > 0) present_[g].push_back(jj);
        }
    }
    a_.assign(N, 0.0); b_.assign(N, 0.0);
    ra_.assign(N, 0); rb_.assign(N, 0);
    ta_.assign(K, 0.0); tb_.assign(K, 0.0);
    tra_.assign(K, 0); trb_.assign(K, 0);
    push_.assign(2 * static_cast<size_t>(N), 0.0);
    dsum_.assign(d_, 0.0);
    xt_.assign(d_, 0.0);
    part_r_.assign(R_, 0);
    part_d_.assign(d_, 0);
}

// the tree as is: every dummy-taxon leaf is a cluster of one taxon
template <class Src>
void Scorer::copy_from(const Src& src, int g) {
    const int off = src.begin(g), end = src.end(g);
    const int base = C_.nodes(), tree = C_.k;
    C_.real_leaf.resize(C_.real_leaf.size() + R_, -1);
    for (int v = off; v < end; ++v) {
        C_.parent.push_back(src.parent(v) < 0 ? -1 : base + (src.parent(v) - off));
        for (int e = src.cbeg(v); e < src.cend(v); ++e) C_.children.push_back(base + (src.child(e) - off));
        C_.child_off.push_back(static_cast<int>(C_.children.size()));
        LeafUnit u{-1, -1, 0.0};
        if (src.leaf(v)) u = src.unit(v);
        if (u.real >= 0) C_.real_leaf[static_cast<size_t>(tree) * R_ + u.real] = base + (v - off);
        C_.leaf_real.push_back(u.real);
        C_.leaf_dummy.push_back(u.dummy);
        C_.leaf_w.push_back(u.w);
    }
    C_.tree_off.push_back(C_.nodes());
    ++C_.k;
}

template <class Src>
void Scorer::compress_from(const Src& src, int g) {
    const int off = src.begin(g), m = src.end(g) - off, root = m - 1;
    const int tag = g + 1;
    auto& down = lab_down_;
    auto& up = lab_up_;
    down.assign(m, kEmpty);
    up.assign(m, kEmpty);
    wsum_.assign(m, 0.0);

    // purity below each node, and units present
    int units = 0;
    for (int v = 0; v < m; ++v) {
        if (src.leaf(off + v)) {
            const LeafUnit u = src.unit(off + v);
            wsum_[v] = u.w;
            if (u.real >= 0) { down[v] = kMixed; ++units; }
            else {
                down[v] = u.dummy;
                if (stamp_[u.dummy] != tag) { stamp_[u.dummy] = tag; ++units; }
            }
        } else {
            int acc = kEmpty;
            double ws = 0;
            for (int e = src.cbeg(off + v); e < src.cend(off + v); ++e) {
                const int c = src.child(e) - off;
                acc = combine(acc, down[c]);
                ws += wsum_[c];
            }
            down[v] = acc;
            wsum_[v] = ws;
        }
    }
    if (units < 4) return;                       // no quartet with 4 distinct units: contributes nothing
    const double W = wsum_[root];

    // purity above each node (the rest of the unrooted tree)
    std::vector<int>& suffix = suffix_;
    for (int v = root; v >= 0; --v) {
        const int cb = src.cbeg(off + v), ce = src.cend(off + v), nc = ce - cb;
        if (nc == 0) continue;
        suffix.assign(nc + 1, kEmpty);
        for (int i = nc - 1; i >= 0; --i) suffix[i] = combine(suffix[i + 1], down[src.child(cb + i) - off]);
        int prefix = up[v];
        for (int i = 0; i < nc; ++i) {
            const int c = src.child(cb + i) - off;
            up[c] = combine(prefix, suffix[i + 1]);
            prefix = combine(prefix, down[c]);
        }
    }

    // kept nodes: not inside a single-dummy component (no dummy has D-1 pure branches)
    const int root_children = src.cend(off + root) - src.cbeg(off + root);
    const bool root_suppressed = root_children == 2;     // a degree-2 root is not a node of the unrooted tree
    kept_.assign(m, 0);
    std::vector<int>& touched = touched_;
    touched.clear();
    int first_kept = -1;
    for (int v = 0; v < m; ++v) {
        const int cb = src.cbeg(off + v), ce = src.cend(off + v);
        if (cb == ce || (v == root && root_suppressed)) continue;
        int D = 0, max_run = 0;
        auto count = [&](int lab) {
            ++D;
            if (lab < 0) return;
            if (slot_[lab] < 0) { slot_[lab] = 0; touched.push_back(lab); }
            max_run = std::max(max_run, ++slot_[lab]);
        };
        for (int e = cb; e < ce; ++e) count(down[src.child(e) - off]);
        if (v != root) count(up[v]);
        for (int j : touched) slot_[j] = -1;
        touched.clear();
        if (D >= 3 && max_run < D - 1) {
            kept_[v] = 1;
            if (first_kept < 0) first_kept = v;
        }
    }
    if (first_kept < 0) { copy_from(src, g); return; }   // not expected for a tree with >= 4 units

    // compressed unrooted tree: kept nodes, real leaves and merged clusters
    int na = 0;
    auto new_node = [&](int r, int dm, double w) {
        if (na < static_cast<int>(adj_.size())) adj_[na].clear(); else adj_.emplace_back();
        if (na < static_cast<int>(node_real_.size())) { node_real_[na] = r; node_dummy_[na] = dm; node_w_[na] = w; }
        else { node_real_.push_back(r); node_dummy_.push_back(dm); node_w_.push_back(w); }
        return na++;
    };
    auto edge = [&](int a, int b) { adj_[a].push_back(b); adj_[b].push_back(a); };
    cid_.assign(m, -1);
    for (int v = 0; v < m; ++v) if (kept_[v]) cid_[v] = new_node(-1, -1, 0);
    auto require_kept = [&](int v) {
        if (!kept_[v]) throw std::logic_error("compression: mixed neighbour of a kept node is not kept");
        return cid_[v];
    };
    std::vector<std::pair<int, double>>& clusters = clusters_;
    clusters.clear();
    for (int v = 0; v < m; ++v) {
        if (!kept_[v]) continue;
        const int u = cid_[v];
        auto add_cluster = [&](int j, double w) {
            if (slot_[j] < 0) { slot_[j] = static_cast<int>(clusters.size()); clusters.emplace_back(j, 0.0); }
            clusters[slot_[j]].second += w;
        };
        for (int e = src.cbeg(off + v); e < src.cend(off + v); ++e) {
            const int c = src.child(e) - off;
            if (down[c] >= 0) add_cluster(down[c], wsum_[c]);
            else if (src.leaf(off + c)) edge(u, new_node(src.unit(off + c).real, -1, 1));
            else edge(u, require_kept(c));
        }
        if (v != root) {
            if (up[v] >= 0) {
                add_cluster(up[v], W - wsum_[v]);
            } else {
                const int p = src.parent(off + v) - off;
                if (p == root && root_suppressed) {          // the edge across the degree-2 root
                    const int cb = src.cbeg(off + root);
                    const int o = src.child(cb) - off == v ? src.child(cb + 1) - off : src.child(cb) - off;
                    if (src.leaf(off + o)) edge(u, new_node(src.unit(off + o).real, -1, 1));
                    else if (v < o) edge(u, require_kept(o));
                }                                             // otherwise the parent (kept) adds the edge
            }
        }
        for (const auto& [j, w] : clusters) { edge(u, new_node(-1, j, w)); slot_[j] = -1; }
        clusters.clear();
    }

    // emit rooted at a kept node, children before parents
    const int tree = C_.k;
    C_.real_leaf.resize(C_.real_leaf.size() + R_, -1);
    ni_.assign(na, -1);
    struct Frame { int x, parent; size_t next; };
    std::vector<Frame> st{{cid_[first_kept], -1, 0}};
    while (!st.empty()) {
        Frame& f = st.back();
        if (f.next < adj_[f.x].size()) {
            const int y = adj_[f.x][f.next++];
            if (y != f.parent) st.push_back({y, f.x, 0});   // invalidates f
            continue;
        }
        const int x = f.x, px = f.parent;
        st.pop_back();
        const int id = C_.nodes();
        ni_[x] = id;
        C_.parent.push_back(-1);
        for (int y : adj_[x]) {
            if (y == px) continue;
            C_.children.push_back(ni_[y]);
            C_.parent[ni_[y]] = id;
        }
        C_.child_off.push_back(static_cast<int>(C_.children.size()));
        C_.leaf_real.push_back(node_real_[x]);
        C_.leaf_dummy.push_back(node_dummy_[x]);
        C_.leaf_w.push_back(node_dummy_[x] >= 0 || node_real_[x] >= 0 ? node_w_[x] : 0.0);
        if (node_real_[x] >= 0) C_.real_leaf[static_cast<size_t>(tree) * R_ + node_real_[x]] = id;
    }
    C_.tree_off.push_back(C_.nodes());
    ++C_.k;
}

void Scorer::set_partition(const std::vector<int8_t>& part_r, const std::vector<int8_t>& part_d) {
    part_r_ = part_r;
    part_d_ = part_d;
    const int N = C_.nodes();
    for (int v = 0; v < N; ++v) {                 // children precede parents
        if (C_.leaf_real[v] >= 0) {
            const int p = part_r_[C_.leaf_real[v]];
            a_[v] = p == 0; b_[v] = p == 1;
            ra_[v] = p == 0; rb_[v] = p == 1;
        } else if (C_.leaf_dummy[v] >= 0) {
            const int p = part_d_[C_.leaf_dummy[v]];
            a_[v] = p == 0 ? C_.leaf_w[v] : 0.0;
            b_[v] = p == 1 ? C_.leaf_w[v] : 0.0;
            ra_[v] = rb_[v] = 0;
        } else {
            double a = 0, b = 0;
            int ra = 0, rb = 0;
            for (int e = C_.child_off[v]; e < C_.child_off[v + 1]; ++e) {
                const int c = C_.children[e];
                a += a_[c]; b += b_[c]; ra += ra_[c]; rb += rb_[c];
            }
            a_[v] = a; b_[v] = b; ra_[v] = ra; rb_[v] = rb;
        }
    }
    for (int g = 0; g < C_.k; ++g) {
        const int root = C_.tree_off[g + 1] - 1;
        ta_[g] = a_[root]; tb_[g] = b_[root];
        tra_[g] = ra_[root]; trb_[g] = rb_[root];
    }
}

double Scorer::score() {
    double s = 0;
    for (int g = 0; g < C_.k; ++g) s += tree_pass(g, false, nullptr, nullptr);
    return s;
}

double Scorer::gains(std::vector<double>& real_gain, std::vector<double>& dummy_gain) {
    real_gain.assign(R_, 0.0);
    dummy_gain.assign(d_, 0.0);
    double s = 0;
    for (int g = 0; g < C_.k; ++g) s += tree_pass(g, true, real_gain.data(), dummy_gain.data());
    return s;
}

// f at the current generic node (branch values in ga_/gb_/grb_/gx_), with dummy `moved` (or none)
// switched to the other side
double Scorer::generic_value(int D, int moved, PolyEval& pe, const std::vector<int>& pres) {
    int mA = 1, mB = 1;
    for (int j : pres) ((part_d_[j] ^ (j == moved)) == 0 ? mA : mB)++;
    pe.resize(D, mA, mB);
    int ia = 1, ib = 1;
    for (int j : pres) {
        const bool in_a = (part_d_[j] ^ (j == moved)) == 0;
        double* dst = in_a ? &pe.al[static_cast<size_t>(ia++) * D] : &pe.be[static_cast<size_t>(ib++) * D];
        for (int i = 0; i < D; ++i) dst[i] = gx_[i][j];
    }
    for (int i = 0; i < D; ++i) {
        double a = ga_[i], b = gb_[i];
        if (moved >= 0) {
            const double x = gx_[i][moved];
            if (part_d_[moved] == 0) { a -= x; b += x; } else { a += x; b -= x; }
        }
        pe.al[i] = a;
        pe.be[i] = b;
        pe.rho[i] = grb_[i];
    }
    pe.build();
    return pe.eval();
}

// `want` is explicit: with no real (or no dummy) taxa the gain vectors are empty and data() may be null
double Scorer::tree_pass(int g, bool want, double* real_gain, double* dummy_gain) {
    const int off = C_.tree_off[g], end = C_.tree_off[g + 1], root = end - 1;
    const double A = ta_[g], B = tb_[g];
    const int RA = tra_[g], RB = trb_[g];
    const double* X = d_ ? &tx_[static_cast<size_t>(g) * d_] : nullptr;
    const std::vector<int>& pres = present_[g];

    double fsum = 0, G0 = 0, G1 = 0;
    if (want) std::fill(dsum_.begin(), dsum_.end(), 0.0);

    for (int v = off; v < end; ++v) {
        const int cb = C_.child_off[v], ce = C_.child_off[v + 1], nc = ce - cb;
        if (nc == 0) continue;
        const double at = A - a_[v], bt = B - b_[v];
        const int rat = RA - ra_[v], rbt = RB - rb_[v];
        const double* xv = d_ ? &x_[static_cast<size_t>(v) * d_] : nullptr;
        for (int j : pres) xt_[j] = X[j] - xv[j];

        if (nc == 2) {
            // branches 0, 1 = children, 2 = top; PA[i] pairs branches (i, i+1), f = sum_i PA[i] PB[i+2]
            const int c0 = C_.children[cb], c1 = C_.children[cb + 1];
            const double a[3] = {a_[c0], a_[c1], at}, b[3] = {b_[c0], b_[c1], bt};
            const int ra3[3] = {ra_[c0], ra_[c1], rat}, rb3[3] = {rb_[c0], rb_[c1], rbt};
            const double* x[3] = {d_ ? &x_[static_cast<size_t>(c0) * d_] : nullptr,
                                  d_ ? &x_[static_cast<size_t>(c1) * d_] : nullptr, xt_.data()};
            double sA[3] = {0, 0, 0}, sB[3] = {double(rb3[0]), double(rb3[1]), double(rb3[2])};
            for (int j : pres) {
                const double x0 = x[0][j], x1 = x[1][j], x2 = x[2][j];
                if (part_d_[j] == 0) { sA[0] += x0 * x1; sA[1] += x1 * x2; sA[2] += x2 * x0; }
                else { sB[0] += x0 * x0; sB[1] += x1 * x1; sB[2] += x2 * x2; }
            }
            double PA[3], PB[3];
            for (int i = 0; i < 3; ++i) {
                PA[i] = a[i] * a[(i + 1) % 3] - sA[i];
                PB[i] = (b[i] * b[i] - sB[i]) / 2;
            }
            const double f = PA[0] * PB[2] + PA[1] * PB[0] + PA[2] * PB[1];
            fsum += f;
            if (!want) continue;

            double dAB[3], dBA[3];
            for (int i = 0; i < 3; ++i) {
                const int j = (i + 1) % 3, l = (i + 2) % 3;
                // the pair opposite branch i is PA[j]; pairs (i,j) and (l,i) lose/gain a_j and a_l
                dAB[i] = ra3[i] > 0 ? PA[j] * b[i] - a[j] * PB[l] - a[l] * PB[j] : 0.0;
                dBA[i] = rb3[i] > 0 ? -PA[j] * (b[i] - 1) + a[j] * PB[l] + a[l] * PB[j] : 0.0;
            }
            push_[2 * c0] = dAB[0] - dAB[2];
            push_[2 * c0 + 1] = dBA[0] - dBA[2];
            push_[2 * c1] = dAB[1] - dAB[2];
            push_[2 * c1 + 1] = dBA[1] - dBA[2];
            G0 += dAB[2];
            G1 += dBA[2];

            for (int jd : pres) {
                const double xs[3] = {x[0][jd], x[1][jd], x[2][jd]};
                const double s = part_d_[jd] == 0 ? -1.0 : 1.0;   // A->B: a -= x, b += x
                double a2[3], b2[3], sA2[3], sB2[3];
                for (int i = 0; i < 3; ++i) {
                    a2[i] = a[i] + s * xs[i];
                    b2[i] = b[i] - s * xs[i];
                    sA2[i] = sA[i] + s * xs[i] * xs[(i + 1) % 3];
                    sB2[i] = sB[i] - s * xs[i] * xs[i];
                }
                double f2 = 0;
                for (int i = 0; i < 3; ++i) {
                    const double pa = a2[i] * a2[(i + 1) % 3] - sA2[i];
                    const double pb = (b2[(i + 2) % 3] * b2[(i + 2) % 3] - sB2[(i + 2) % 3]) / 2;
                    f2 += pa * pb;
                }
                dsum_[jd] += f2 - f;
            }
        } else if (pres.empty()) {
            // polytomy in a tree without dummy taxa: fixed-size power sums
            Agg0 g{};
            for (int k = 0; k < nc; ++k) agg0_add(g, a_[C_.children[cb + k]], b_[C_.children[cb + k]], 1);
            agg0_add(g, at, bt, 1);
            const double f0 = agg0_f(g);
            fsum += f0;
            if (!want) continue;
            auto delta = [&](double a, double b, int dir) {
                Agg0 h = g;
                agg0_add(h, a, b, -1);
                if (dir == 0) agg0_add(h, a - 1, b + 1, 1); else agg0_add(h, a + 1, b - 1, 1);
                return agg0_f(h) - f0;
            };
            const double tAB = rat > 0 ? delta(at, bt, 0) : 0.0, tBA = rbt > 0 ? delta(at, bt, 1) : 0.0;
            for (int k = 0; k < nc; ++k) {
                const int c = C_.children[cb + k];
                push_[2 * c] = (ra_[c] > 0 ? delta(a_[c], b_[c], 0) : 0.0) - tAB;
                push_[2 * c + 1] = (rb_[c] > 0 ? delta(a_[c], b_[c], 1) : 0.0) - tBA;
            }
            G0 += tAB;
            G1 += tBA;
        } else {
            const int D = nc + 1;
            ga_.resize(D); gb_.resize(D); gra_.resize(D); grb_.resize(D); gx_.resize(D); gd_.resize(2 * D);
            for (int k = 0; k < nc; ++k) {
                const int c = C_.children[cb + k];
                ga_[k] = a_[c]; gb_[k] = b_[c]; gra_[k] = ra_[c]; grb_[k] = rb_[c];
                gx_[k] = d_ ? &x_[static_cast<size_t>(c) * d_] : nullptr;
            }
            ga_[nc] = at; gb_[nc] = bt; gra_[nc] = rat; grb_[nc] = rbt; gx_[nc] = xt_.data();

            const double f0 = generic_value(D, -1, pe_, pres);
            fsum += f0;
            if (!want) continue;

            agg_backup_ = pe_.agg;
            for (int i = 0; i < D; ++i) {
                for (int dir = 0; dir < 2; ++dir) {
                    double delta = 0;
                    if ((dir == 0 ? gra_[i] : grb_[i]) > 0) {
                        const double al0 = pe_.al[i], be0 = pe_.be[i], rho0 = pe_.rho[i];
                        const double s = dir == 0 ? -1.0 : 1.0;   // A->B: a-1, b+1, rho+1
                        pe_.add(i, -1);
                        pe_.al[i] = al0 + s;
                        pe_.be[i] = be0 - s;
                        pe_.rho[i] = rho0 - s;
                        pe_.add(i, 1);
                        delta = pe_.eval() - f0;
                        pe_.al[i] = al0; pe_.be[i] = be0; pe_.rho[i] = rho0;
                        pe_.agg = agg_backup_;
                    }
                    gd_[2 * i + dir] = delta;
                }
            }
            for (int k = 0; k < nc; ++k) {
                const int c = C_.children[cb + k];
                push_[2 * c] = gd_[2 * k] - gd_[2 * nc];
                push_[2 * c + 1] = gd_[2 * k + 1] - gd_[2 * nc + 1];
            }
            G0 += gd_[2 * nc];
            G1 += gd_[2 * nc + 1];
            for (int jd : pres) dsum_[jd] += generic_value(D, jd, pe2_, pres) - f0;
        }
    }

    // per-tree totals: w(PA(g)) w(PB(g)) = all quartets with two taxa on each side (Section 2.5.2)
    double sqA = 0, sqB = 0;
    for (int j : pres) (part_d_[j] == 0 ? sqA : sqB) += X[j] * X[j];
    const double PA = (A * A - RA - sqA) / 2, PB = (B * B - RB - sqB) / 2, tq = PA * PB;
    const double score = 2 * fsum - tq;
    if (!want) return score;

    // real taxa: gain = path sum of pushed deltas + top-branch deltas of every node
    push_[2 * root] = 0;
    push_[2 * root + 1] = 0;
    const double dtot0 = (PA - (A - 1)) * (PB + B) - tq;   // one real taxon A -> B
    const double dtot1 = (PA + A) * (PB - (B - 1)) - tq;   // B -> A
    for (int v = root - 1; v >= off; --v) {
        const int p = C_.parent[v];
        push_[2 * v] += push_[2 * p];
        push_[2 * v + 1] += push_[2 * p + 1];
        const int i = C_.leaf_real[v];
        if (i >= 0)
            real_gain[i] += part_r_[i] == 0 ? 2 * (push_[2 * v] + G0) - dtot0 : 2 * (push_[2 * v + 1] + G1) - dtot1;
    }
    for (int j : pres) {
        const double x = X[j];
        const double dt = part_d_[j] == 0 ? (PA - x * (A - x)) * (PB + x * B) - tq
                                          : (PA + x * A) * (PB - x * (B - x)) - tq;
        dummy_gain[j] += 2 * dsum_[j] - dt;
    }
    return score;
}

void Scorer::swap_real(int i) {
    const int p = part_r_[i];
    for (int g = 0; g < C_.k; ++g) {
        int v = C_.real_leaf[static_cast<size_t>(g) * R_ + i];
        if (v < 0) continue;
        for (; v >= 0; v = C_.parent[v]) {
            if (p == 0) { a_[v] -= 1; b_[v] += 1; --ra_[v]; ++rb_[v]; }
            else { a_[v] += 1; b_[v] -= 1; ++ra_[v]; --rb_[v]; }
        }
        if (p == 0) { ta_[g] -= 1; tb_[g] += 1; --tra_[g]; ++trb_[g]; }
        else { ta_[g] += 1; tb_[g] -= 1; ++tra_[g]; --trb_[g]; }
    }
    part_r_[i] = static_cast<int8_t>(1 - p);
}

void Scorer::swap_dummy(int j) {
    const int p = part_d_[j];
    const int N = C_.nodes();
    for (int v = 0; v < N; ++v) {
        const double x = x_[static_cast<size_t>(v) * d_ + j];
        if (x == 0) continue;
        if (p == 0) { a_[v] -= x; b_[v] += x; } else { a_[v] += x; b_[v] -= x; }
    }
    for (int g = 0; g < C_.k; ++g) {
        const double x = tx_[static_cast<size_t>(g) * d_ + j];
        if (p == 0) { ta_[g] -= x; tb_[g] += x; } else { ta_[g] += x; tb_[g] -= x; }
    }
    part_d_[j] = static_cast<int8_t>(1 - p);
}

}  // namespace wqfm
