#include "scorer.hpp"

#include <algorithm>
#include <stdexcept>

namespace wqfm {

// ---------------------------------------------------------------- PolyNode
//
// Components: m = 0 is all of A (alpha_0i = a_i) or all of B (beta_0i = b_i), with sign +1; each dummy
// is a component of its side with sign -1. With per-component sums over the branches TA = sum alpha,
// QA = sum alpha^2, R = sum pb alpha, R2 = sum pb alpha^2 (pb_i = PB_i), TB, QB likewise, P = sum pb, and
// per pair (m in A, n in B) M1..M4 = sum alpha beta, alpha^2 beta, alpha beta^2, alpha^2 beta^2:
//   EA = sum_m s_m (TA_m^2 - QA_m)/2,  EB likewise,  S = P EA - sum_m s_m (TA_m R_m - R2_m),
//   U  = EA EB - sum_{m,n} s_m s_n tau(TA_m, TB_n, M_mn),
//   tau(ta, tb, M) = tb ta M1 - tb M2 - ta M3 + M4 - (M1^2 - M4)/2.
// Pairs of components that share no branch have M = 0 and contribute nothing. The pair (A, B) uses N1..N4,
// (A, dummy l in B) uses G_l = sums of a x_l, a^2 x_l, a x_l^2, a^2 x_l^2, (dummy k in A, B) uses H_k =
// sums of x_k b, x_k^2 b, x_k b^2, x_k^2 b^2, and their sums over each side are linear in TA (TB):
//   sum_{l in B} tau(TA, T_l, G_l) = TA CG1 + CG0,  sum_{k in A} tau(T_k, TB, H_k) = TB CH1 + CH0.
// Dummy pairs depend on the dummies' sides only; their total over opposite-side pairs is Xp.

namespace {
inline double tau(double ta, double tb, double m1, double m2, double m3, double m4) {
    return tb * ta * m1 - tb * m2 - ta * m3 + m4 - (m1 * m1 - m4) / 2;
}
}  // namespace

void PolyNode::begin(int dummies) {
    p_ = dummies;
    a_.clear(); b_.clear(); rho_.clear(); sq0_.clear();
    eoff_.assign(1, 0);
    ek_.clear(); ex_.clear();
    inA.resize(dummies);
}

void PolyNode::branch(double a, double b, double rho, double sq_extra) {
    if (!a_.empty()) eoff_.push_back(static_cast<int>(ek_.size()));
    a_.push_back(a);
    b_.push_back(b);
    rho_.push_back(rho);
    sq0_.push_back(sq_extra);
}

void PolyNode::entry(int k, double x) {
    ek_.push_back(k);
    ex_.push_back(x);
}

double PolyNode::eval(const Sums& s) {
    const double EA = (s.TA0 * s.TA0 - s.QA0) / 2 - s.EAd, EB = (s.TB0 * s.TB0 - s.QB0) / 2 - s.EBd;
    const double sr = (s.TA0 * s.R0 - s.R20) - s.srd;
    const double X = tau(s.TA0, s.TB0, s.N1, s.N2, s.N3, s.N4) - (s.TA0 * s.CG1 + s.CG0) - (s.TB0 * s.CH1 + s.CH0) + s.Xp;
    return s.P * EA - sr + (EA * EB - X) / 2;
}

// add (sign +1) or remove (-1) the side-total terms of a dummy with sums t, q and g (R R2 G1..G4 H1..H4)
void PolyNode::terms(Sums& s, bool in_a, double t, double q, const double* g, double sign) {
    if (in_a) {
        s.EAd += sign * (t * t - q) / 2;
        s.srd += sign * (t * g[0] - g[1]);
        s.CH1 += sign * (t * g[6] - g[7]);
        s.CH0 += sign * (-t * g[8] + g[9] - (g[6] * g[6] - g[9]) / 2);
    } else {
        s.EBd += sign * (t * t - q) / 2;
        s.CG1 += sign * (t * g[2] - g[4]);
        s.CG0 += sign * (-t * g[3] + g[5] - (g[2] * g[2] - g[5]) / 2);
    }
}

double PolyNode::build() {
    const int D = static_cast<int>(a_.size()), p = p_;
    eoff_.push_back(static_cast<int>(ek_.size()));
    sqB_.assign(sq0_.begin(), sq0_.end());
    pb_.resize(D);
    T_.assign(p, 0.0);
    Q_.assign(p, 0.0);
    agg_.assign(static_cast<size_t>(p) * kAgg, 0.0);
    koff_.assign(p + 1, 0);
    s_ = Sums{};
    for (int i = 0; i < D; ++i) {
        for (int e = eoff_[i]; e < eoff_[i + 1]; ++e) {
            if (!inA[ek_[e]]) sqB_[i] += ex_[e] * ex_[e];
            ++koff_[ek_[e] + 1];
        }
        const double a = a_[i], b = b_[i], pb = (b * b - sqB_[i] - rho_[i]) / 2;
        pb_[i] = pb;
        s_.TA0 += a; s_.QA0 += a * a; s_.R0 += pb * a; s_.R20 += pb * a * a;
        s_.TB0 += b; s_.QB0 += b * b; s_.P += pb;
        s_.N1 += a * b; s_.N2 += a * a * b; s_.N3 += a * b * b; s_.N4 += a * a * b * b;
        for (int e = eoff_[i]; e < eoff_[i + 1]; ++e) {
            const int k = ek_[e];
            const double x = ex_[e];
            double* g = &agg_[static_cast<size_t>(k) * kAgg];
            T_[k] += x; Q_[k] += x * x;
            g[0] += pb * x; g[1] += pb * x * x;
            g[2] += a * x; g[3] += a * a * x; g[4] += a * x * x; g[5] += a * a * x * x;
            g[6] += x * b; g[7] += x * x * b; g[8] += x * b * b; g[9] += x * x * b * b;
        }
    }
    // each dummy's branches
    for (int k = 0; k < p; ++k) koff_[k + 1] += koff_[k];
    kbr_.resize(koff_[p]);
    kx_.resize(koff_[p]);
    kfill_.assign(koff_.begin(), koff_.end() - 1);
    for (int i = 0; i < D; ++i)
        for (int e = eoff_[i]; e < eoff_[i + 1]; ++e) {
            const int at = kfill_[ek_[e]]++;
            kbr_[at] = i;
            kx_[at] = ex_[e];
        }
    // dummy terms per side, and each dummy's pair terms against the other side as a member of A or of B
    for (int k = 0; k < p; ++k) terms(s_, inA[k] != 0, T_[k], Q_[k], &agg_[static_cast<size_t>(k) * kAgg], 1);
    pair_terms(D);
    for (int k = 0; k < p; ++k) if (inA[k]) s_.Xp += pairA_[k];
    if (dacc_.size() < static_cast<size_t>(p) * kAgg) dacc_.resize(static_cast<size_t>(p) * kAgg);   // zeroed on first use
    dstamp_.assign(p, 0);
    tag_ = 0;
    f0_ = eval(s_);
    return f0_;
}

// Pair terms of the dummies from sums over each side, in O(entries) instead of O(p^2) over pairs (at a
// node nearly every pair of listed dummies shares the top branch). With m_kl the sums over the branches
// of x_k x_l, x_k^2 x_l, x_k x_l^2 and x_k^2 x_l^2, the term of k in A and l in B is
//   tau(T_k, T_l, m_kl) = T_k T_l m1 - T_l m2 - T_k m3 + 3/2 m4 - m1^2 / 2.
// Summed over l in B, each term but m1^2 is a sum over k's branches i of x_ki or x_ki^2 times
// B1_i = sum_l T_l x_li or B2_i = sum_l x_li^2, and sum_l m1_kl^2 = sum_ij x_ki x_kj G_ij with
// G_ij = sum_l x_li x_lj: its diagonal is B2_i, its other entries are kept for the pairs of branches one
// dummy spans. So pairA_[k] = T_k sB1 - sB1q - T_k sB2 + sB2q - offB, where s.. sum over k's branches
// (q: with x_ki^2) and offB = sum_{i<j} x_ki x_kj G_ij; pairB_[l] likewise over A, with the roles of
// m2 and m3 swapped. A dummy on the summed side is in these sums itself, so its term with itself,
// T^2 M1 - 2 T M2 + 3/2 M4 - M1^2 / 2 (M_a = sum_i x_ki^a), is subtracted.
void PolyNode::pair_terms(int D) {
    const int p = p_;
    // per dummy, over its branches i (ascending): x_ki and x_ki^2 times the branch's side totals, and
    // the power sums M1 M2 M4 (kPair each); the side totals of a branch come from its entries, so the
    // cost is O(entries), not O(D) as with per-branch arrays (polytomies can have hundreds of branches)
    acc_.assign(static_cast<size_t>(p) * kPair, 0.0);
    for (int i = 0; i < D; ++i) {
        const int eb = eoff_[i], ee = eoff_[i + 1];
        if (eb == ee) continue;
        double A1 = 0, A2 = 0, B1 = 0, B2 = 0;
        for (int e = eb; e < ee; ++e) {
            const int k = ek_[e];
            const double x = ex_[e];
            if (inA[k]) { A1 += T_[k] * x; A2 += x * x; }
            else { B1 += T_[k] * x; B2 += x * x; }
        }
        for (int e = eb; e < ee; ++e) {
            double* s = &acc_[static_cast<size_t>(ek_[e]) * kPair];
            const double x = ex_[e], x2 = x * x;
            s[0] += x * A1; s[1] += x2 * A1; s[2] += x * A2; s[3] += x2 * A2;
            s[4] += x * B1; s[5] += x2 * B1; s[6] += x * B2; s[7] += x2 * B2;
            s[8] += x2; s[9] += x2 * x; s[10] += x2 * x2;
        }
    }
    // off-diagonal G per side: entries (i*D + j, dummy) summed in that order
    gent_.clear();
    for (int k = 0; k < p; ++k)
        for (int u = koff_[k]; u < koff_[k + 1]; ++u)
            for (int w = u + 1; w < koff_[k + 1]; ++w)
                gent_.push_back({static_cast<std::int64_t>(kbr_[u]) * D + kbr_[w], k, kx_[u] * kx_[w]});
    std::sort(gent_.begin(), gent_.end(),
              [](const GramEntry& x, const GramEntry& y) { return x.key != y.key ? x.key < y.key : x.k < y.k; });
    gkey_.clear();
    gA_.clear();
    gB_.clear();
    for (const GramEntry& g : gent_) {
        if (gkey_.empty() || gkey_.back() != g.key) {
            gkey_.push_back(g.key);
            gA_.push_back(0.0);
            gB_.push_back(0.0);
        }
        (inA[g.k] ? gA_.back() : gB_.back()) += g.v;
    }
    pairA_.resize(p);
    pairB_.resize(p);
    for (int k = 0; k < p; ++k) {
        const double* s = &acc_[static_cast<size_t>(k) * kPair];
        const double sA1 = s[0], sA1q = s[1], sA2 = s[2], sA2q = s[3], sB1 = s[4], sB1q = s[5], sB2 = s[6],
                     sB2q = s[7], M1 = s[8], M2 = s[9], M4 = s[10];
        double offA = 0, offB = 0;
        for (int u = koff_[k]; u < koff_[k + 1]; ++u)
            for (int w = u + 1; w < koff_[k + 1]; ++w) {
                const std::int64_t key = static_cast<std::int64_t>(kbr_[u]) * D + kbr_[w];
                const size_t s = static_cast<size_t>(std::lower_bound(gkey_.begin(), gkey_.end(), key) - gkey_.begin());
                const double xx = kx_[u] * kx_[w];
                offA += xx * gA_[s];
                offB += xx * gB_[s];
            }
        const double T = T_[k], self = T * T * M1 - 2 * T * M2 + 1.5 * M4 - 0.5 * M1 * M1;
        pairA_[k] = T * sB1 - sB1q - T * sB2 + sB2q - offB - (inA[k] ? 0.0 : self);   // against l in B
        pairB_[k] = T * sA1 - T * sA2 - sA1q + sA2q - offA - (inA[k] ? self : 0.0);   // against k' in A
    }
}

// sums over the branches when branch i takes new values (dummy sums not included)
void PolyNode::branch_change(Sums& s, int i, double a2, double b2, double rho2, double sqB2) const {
    const double a1 = a_[i], b1 = b_[i], p1 = pb_[i];
    const double p2 = (b2 * b2 - sqB2 - rho2) / 2;
    s.TA0 += a2 - a1; s.QA0 += a2 * a2 - a1 * a1;
    s.R0 += p2 * a2 - p1 * a1; s.R20 += p2 * a2 * a2 - p1 * a1 * a1;
    s.TB0 += b2 - b1; s.QB0 += b2 * b2 - b1 * b1; s.P += p2 - p1;
    s.N1 += a2 * b2 - a1 * b1; s.N2 += a2 * a2 * b2 - a1 * a1 * b1;
    s.N3 += a2 * b2 * b2 - a1 * b1 * b1; s.N4 += a2 * a2 * b2 * b2 - a1 * a1 * b1 * b1;
}

namespace {
// change of a dummy's sums (R R2 G1..G4 H1..H4) from its weight x in a branch whose pb, a, b change
inline void sums_change(double* d, double x, double dp, double da, double da2, double db, double db2) {
    d[0] += dp * x; d[1] += dp * x * x;
    d[2] += da * x; d[3] += da2 * x; d[4] += da * x * x; d[5] += da2 * x * x;
    d[6] += x * db; d[7] += x * x * db; d[8] += x * db2; d[9] += x * x * db2;
}
}  // namespace

double PolyNode::pairs_without(int i) const {
    const double EA = (s_.TA0 * s_.TA0 - s_.QA0) / 2 - s_.EAd;
    double with_i = a_[i] * (s_.TA0 - a_[i]);           // sum over j != i of PA_ij
    for (int e = eoff_[i]; e < eoff_[i + 1]; ++e)
        if (inA[ek_[e]]) with_i -= ex_[e] * (T_[ek_[e]] - ex_[e]);
    return EA - with_i;
}

double PolyNode::real_delta(int i, int dir) {
    const double s = dir == 0 ? -1.0 : 1.0;          // A->B: a-1, b+1, rho+1
    const double a2 = a_[i] + s, b2 = b_[i] - s, rho2 = rho_[i] - s;
    Sums t = s_;
    branch_change(t, i, a2, b2, rho2, sqB_[i]);
    const double dp = (b2 * b2 - sqB_[i] - rho2) / 2 - pb_[i];
    const double da = a2 - a_[i], da2 = a2 * a2 - a_[i] * a_[i], db = b2 - b_[i], db2 = b2 * b2 - b_[i] * b_[i];
    for (int e = eoff_[i]; e < eoff_[i + 1]; ++e) {
        const int k = ek_[e];
        const double* g = &agg_[static_cast<size_t>(k) * kAgg];
        double g2[kAgg];
        std::copy(g, g + kAgg, g2);
        sums_change(g2, ex_[e], dp, da, da2, db, db2);
        terms(t, inA[k] != 0, T_[k], Q_[k], g, -1);
        terms(t, inA[k] != 0, T_[k], Q_[k], g2, 1);
    }
    return eval(t) - f0_;
}

double PolyNode::dummy_delta(int k) {
    const bool was_a = inA[k] != 0;
    Sums t = s_;
    terms(t, was_a, T_[k], Q_[k], &agg_[static_cast<size_t>(k) * kAgg], -1);
    t.Xp += was_a ? pairB_[k] - pairA_[k] : pairA_[k] - pairB_[k];
    // its branches change: A->B: a - x, b + x, x^2 joins the B squares (B->A the reverse)
    ++tag_;
    dtouched_.clear();
    for (int u = koff_[k]; u < koff_[k + 1]; ++u) {
        const int i = kbr_[u];
        const double x = was_a ? kx_[u] : -kx_[u];
        const double a2 = a_[i] - x, b2 = b_[i] + x, sqB2 = sqB_[i] + (was_a ? 1 : -1) * kx_[u] * kx_[u];
        branch_change(t, i, a2, b2, rho_[i], sqB2);
        const double dp = (b2 * b2 - sqB2 - rho_[i]) / 2 - pb_[i];
        const double da = a2 - a_[i], da2 = a2 * a2 - a_[i] * a_[i], db = b2 - b_[i], db2 = b2 * b2 - b_[i] * b_[i];
        for (int e = eoff_[i]; e < eoff_[i + 1]; ++e) {
            const int l = ek_[e];
            double* d = &dacc_[static_cast<size_t>(l) * kAgg];
            if (dstamp_[l] != tag_) {
                dstamp_[l] = tag_;
                dtouched_.push_back(l);
                std::fill(d, d + kAgg, 0.0);
            }
            sums_change(d, ex_[e], dp, da, da2, db, db2);
        }
    }
    for (int l : dtouched_) {
        const double* g = &agg_[static_cast<size_t>(l) * kAgg];
        const double* d = &dacc_[static_cast<size_t>(l) * kAgg];
        double g2[kAgg];
        for (int j = 0; j < kAgg; ++j) g2[j] = g[j] + d[j];
        if (l == k) {
            terms(t, !was_a, T_[k], Q_[k], g2, 1);
        } else {
            terms(t, inA[l] != 0, T_[l], Q_[l], g, -1);
            terms(t, inA[l] != 0, T_[l], Q_[l], g2, 1);
        }
    }
    if (dstamp_[k] != tag_) terms(t, !was_a, T_[k], Q_[k], &agg_[static_cast<size_t>(k) * kAgg], 1);
    return eval(t) - f0_;
}

// ---------------------------------------------------------------- Scorer

namespace {

// Node contribution with no dummy taxa in the tree (all weight real, so rho = b): fixed-size power
// sums, the PolyNode formulas with one component on each side.
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
    moves_.clear();
    C_ = SubForest{};
    C_.tree_off.push_back(0);
    C_.child_off.push_back(0);
    stamp_.assign(d_, 0);
    slot_.assign(d_, -1);
}

void Scorer::finish_setup() {
    const int N = C_.nodes(), K = C_.k;
    // per node, the dummies with weight below it: a cluster leaf's own, or the union of its children's
    // (summed in child order)
    xoff_.assign(N + 1, 0);
    xj_.clear();
    xw_.clear();
    xc0_.clear();
    xc1_.clear();
    if (d_ > 0) {
        std::vector<double> acc(d_, 0.0);
        std::vector<int> seen;
        for (int v = 0; v < N; ++v) {
            if (C_.leaf_dummy[v] >= 0) {
                xj_.push_back(C_.leaf_dummy[v]);
                xw_.push_back(C_.leaf_w[v]);
            } else {
                seen.clear();
                for (int e = C_.child_off[v]; e < C_.child_off[v + 1]; ++e) {
                    const int c = C_.children[e];
                    for (int t = xoff_[c]; t < xoff_[c + 1]; ++t) {
                        if (acc[xj_[t]] == 0) seen.push_back(xj_[t]);
                        acc[xj_[t]] += xw_[t];
                    }
                }
                std::sort(seen.begin(), seen.end());
                for (int j : seen) {
                    xj_.push_back(j);
                    xw_.push_back(acc[j]);
                    acc[j] = 0;
                }
            }
            xoff_[v + 1] = static_cast<int>(xj_.size());
        }
        // at binary nodes, each listed dummy's weights below the two children
        xc0_.assign(xj_.size(), 0.0);
        xc1_.assign(xj_.size(), 0.0);
        for (int v = 0; v < N; ++v) {
            if (C_.child_off[v + 1] - C_.child_off[v] != 2) continue;
            for (int side = 0; side < 2; ++side) {
                const int c = C_.children[C_.child_off[v] + side];
                std::vector<double>& out = side == 0 ? xc0_ : xc1_;
                int e = xoff_[v];
                for (int t = xoff_[c]; t < xoff_[c + 1]; ++t) {
                    while (xj_[e] != xj_[t]) ++e;          // the child's list is a sublist of v's
                    out[e] = xw_[t];
                }
            }
        }
    }
    // per dummy, its nodes
    doff_.assign(d_ + 1, 0);
    for (int j : xj_) ++doff_[j + 1];
    for (int j = 0; j < d_; ++j) doff_[j + 1] += doff_[j];
    dnode_.resize(xj_.size());
    dw_.resize(xj_.size());
    {
        std::vector<int> fill(doff_.begin(), doff_.end() - 1);
        for (int v = 0; v < N; ++v)
            for (int t = xoff_[v]; t < xoff_[v + 1]; ++t) {
                const int at = fill[xj_[t]]++;
                dnode_[at] = v;
                dw_[at] = xw_[t];
            }
    }
    tx_.assign(static_cast<size_t>(K) * d_, 0.0);
    present_.assign(K, {});
    scale_ = 0;
    for (int g = 0; g < K; ++g) {
        double w = 0;
        for (int v = C_.tree_off[g]; v < C_.tree_off[g + 1]; ++v)
            if (C_.leaf_real[v] >= 0 || C_.leaf_dummy[v] >= 0) w += C_.leaf_w[v];
        scale_ += w * w * w * w;
        const int root = C_.tree_off[g + 1] - 1;
        for (int t = xoff_[root]; t < xoff_[root + 1]; ++t) {
            tx_[static_cast<size_t>(g) * d_ + xj_[t]] = xw_[t];
            present_[g].push_back(xj_[t]);
        }
    }
    inner_.clear(); rleaf_.clear();
    inner_off_.assign(1, 0); rleaf_off_.assign(1, 0);
    for (int g = 0; g < K; ++g) {
        for (int v = C_.tree_off[g]; v < C_.tree_off[g + 1]; ++v) {
            if (C_.child_off[v] < C_.child_off[v + 1]) inner_.push_back(v);
            else if (C_.leaf_real[v] >= 0) rleaf_.push_back(v);
        }
        inner_off_.push_back(static_cast<int>(inner_.size()));
        rleaf_off_.push_back(static_cast<int>(rleaf_.size()));
    }
    a_.assign(N, 0.0); b_.assign(N, 0.0);
    ra_.assign(N, 0); rb_.assign(N, 0);
    ta_.assign(K, 0.0); tb_.assign(K, 0.0);
    tra_.assign(K, 0); trb_.assign(K, 0);
    push_.assign(2 * static_cast<size_t>(N), 0.0);

    // chunks of whole trees, balanced by node count; their number depends on the forest only
    constexpr int kChunkNodes = 4096, kMaxChunks = 128;
    const int C = std::max(1, std::min({kMaxChunks, K, N / kChunkNodes}));
    chunk_.assign(1, 0);
    for (int c = 1; c < C; ++c) {
        const int target = static_cast<int>(static_cast<long long>(N) * c / C);
        int g = static_cast<int>(std::lower_bound(C_.tree_off.begin(), C_.tree_off.end(), target) - C_.tree_off.begin());
        g = std::max(g, chunk_.back());
        if (g > chunk_.back() && g < K) chunk_.push_back(g);
    }
    chunk_.push_back(K);
    chunk_score_.assign(chunk_.size() - 1, 0.0);
    chunk_gain_.assign((chunk_.size() - 1) * static_cast<size_t>(R_ + d_), 0.0);
    ensure_scratch();
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
    moves_.clear();                               // the totals are rebuilt from part_r_ / part_d_
    run_chunks([&](int c, int) { set_partition_chunk(c, part_r_.data(), part_d_.data()); });
}

void Scorer::set_partition_chunk(int c, const int8_t* pr, const int8_t* pd) {
    for (int g = chunk_[c]; g < chunk_[c + 1]; ++g) {
        const int end = C_.tree_off[g + 1];
        for (int v = C_.tree_off[g]; v < end; ++v) {   // children precede parents
            if (C_.leaf_real[v] >= 0) {
                const int p = pr[C_.leaf_real[v]];
                a_[v] = p == 0; b_[v] = p == 1;
                ra_[v] = p == 0; rb_[v] = p == 1;
            } else if (C_.leaf_dummy[v] >= 0) {
                const int p = pd[C_.leaf_dummy[v]];
                a_[v] = p == 0 ? C_.leaf_w[v] : 0.0;
                b_[v] = p == 1 ? C_.leaf_w[v] : 0.0;
                ra_[v] = rb_[v] = 0;
            } else {
                double a = 0, b = 0;
                int ra = 0, rb = 0;
                for (int e = C_.child_off[v]; e < C_.child_off[v + 1]; ++e) {
                    const int ch = C_.children[e];
                    a += a_[ch]; b += b_[ch]; ra += ra_[ch]; rb += rb_[ch];
                }
                a_[v] = a; b_[v] = b; ra_[v] = ra; rb_[v] = rb;
            }
        }
        const int root = end - 1;
        ta_[g] = a_[root]; tb_[g] = b_[root];
        tra_[g] = ra_[root]; trb_[g] = rb_[root];
    }
}

void Scorer::set_pool(ThreadPool* pool) {
    pool_ = pool;
    ensure_scratch();
}

void Scorer::ensure_scratch() {
    const size_t n = pool_ ? static_cast<size_t>(pool_->threads()) : 1;
    if (scratch_.size() < n) scratch_.resize(n);
    for (Scratch& sc : scratch_) {
        sc.dsum.assign(d_, 0.0);
        sc.pos.assign(d_, 0);
    }
}

void Scorer::run_chunks(const std::function<void(int, int)>& fn) {
    const int C = static_cast<int>(chunk_.size()) - 1;
    if (pool_ && C > 1) pool_->parallel_for(C, fn);
    else for (int c = 0; c < C; ++c) fn(c, 0);
}

double Scorer::score() {
    run_chunks([&](int c, int slot) {
        apply_moves(c);
        double s = 0;
        for (int g = chunk_[c]; g < chunk_[c + 1]; ++g)
            s += tree_pass(g, false, nullptr, nullptr, scratch_[slot], part_r_.data(), part_d_.data());
        chunk_score_[c] = s;
    });
    moves_.clear();
    double s = 0;
    for (double x : chunk_score_) s += x;
    return s;
}

void Scorer::scores(const Batch& batch, std::vector<double>& out) {
    const int Q = batch.size(), C = static_cast<int>(chunk_.size()) - 1, W = R_ + d_;
    out.assign(Q, 0.0);
    if (Q == 0) return;
    batch_score_.resize(static_cast<size_t>(Q) * C);
    run_chunks([&](int c, int slot) {
        apply_moves(c);
        for (int q = 0; q < Q; ++q) {
            const int8_t* pr = batch.sides.data() + static_cast<size_t>(q) * W;
            const int8_t* pd = pr + R_;
            if (batch.reset[q]) set_partition_chunk(c, pr, pd);
            else for (int e = batch.off[q]; e < batch.off[q + 1]; ++e) move_real(c, batch.moved[e], 1 - pr[batch.moved[e]]);
            double s = 0;
            for (int g = chunk_[c]; g < chunk_[c + 1]; ++g) s += tree_pass(g, false, nullptr, nullptr, scratch_[slot], pr, pd);
            batch_score_[static_cast<size_t>(q) * C + c] = s;
        }
    });
    moves_.clear();
    const int8_t* last = batch.sides.data() + static_cast<size_t>(Q - 1) * W;
    part_r_.assign(last, last + R_);
    part_d_.assign(last + R_, last + W);
    for (int q = 0; q < Q; ++q) {               // chunk order, as in score()
        double s = 0;
        for (int c = 0; c < C; ++c) s += batch_score_[static_cast<size_t>(q) * C + c];
        out[q] = s;
    }
}

double Scorer::gains(std::vector<double>& real_gain, std::vector<double>& dummy_gain) {
    const int W = R_ + d_;
    run_chunks([&](int c, int slot) {
        apply_moves(c);
        double* rg = &chunk_gain_[static_cast<size_t>(c) * W];
        std::fill(rg, rg + W, 0.0);
        double s = 0;
        for (int g = chunk_[c]; g < chunk_[c + 1]; ++g)
            s += tree_pass(g, true, rg, rg + R_, scratch_[slot], part_r_.data(), part_d_.data());
        chunk_score_[c] = s;
    });
    moves_.clear();
    real_gain.assign(R_, 0.0);
    dummy_gain.assign(d_, 0.0);
    double s = 0;
    for (size_t c = 0; c < chunk_score_.size(); ++c) {    // fixed order: same result for any thread count
        s += chunk_score_[c];
        const double* rg = &chunk_gain_[c * W];
        for (int i = 0; i < R_; ++i) real_gain[i] += rg[i];
        for (int j = 0; j < d_; ++j) dummy_gain[j] += rg[R_ + j];
    }
    return s;
}

// tree_pass for a tree with no dummy-taxon leaves (every tree at the top level). All leaves are real
// taxa of weight 1, so a and b equal the real-taxon counts ra and rb, the dummy terms vanish, and
// polytomies use fixed-size power sums (Agg0). The binary-node arithmetic is that of
// tree_pass_dummies with the dummy terms dropped.
double Scorer::tree_pass_plain(int g, bool want, double* real_gain, const int8_t* pr) {
    const int root = C_.tree_off[g + 1] - 1;
    const double A = ta_[g], B = tb_[g];
    const int RA = tra_[g], RB = trb_[g];
    const int* __restrict co = C_.child_off.data();
    const int* __restrict kids = C_.children.data();
    const double* __restrict av = a_.data();
    const double* __restrict bv = b_.data();
    double* __restrict push = push_.data();

    const int* __restrict inner = inner_.data();
    const int qb = inner_off_[g], qe = inner_off_[g + 1];

    double fsum = 0, G0 = 0, G1 = 0;
    for (int q = qb; q < qe; ++q) {
        const int v = inner[q];
        const int cb = co[v], nc = co[v + 1] - cb;
        const double at = A - av[v], bt = B - bv[v];
        if (nc == 2) {
            const int c0 = kids[cb], c1 = kids[cb + 1];
            const double a0 = av[c0], a1 = av[c1], b0 = bv[c0], b1 = bv[c1];
            const double PA0 = a0 * a1, PA1 = a1 * at, PA2 = at * a0;
            const double PB0 = (b0 * b0 - b0) / 2, PB1 = (b1 * b1 - b1) / 2, PB2 = (bt * bt - bt) / 2;
            fsum += PA0 * PB2 + PA1 * PB0 + PA2 * PB1;
            if (!want) continue;
            // dAB[i] = PA[i+1] b[i] - a[i+1] PB[i+2] - a[i+2] PB[i+1], dBA likewise (see tree_pass)
            const double dAB0 = a0 > 0 ? PA1 * b0 - a1 * PB2 - at * PB1 : 0.0;
            const double dAB1 = a1 > 0 ? PA2 * b1 - at * PB0 - a0 * PB2 : 0.0;
            const double dAB2 = at > 0 ? PA0 * bt - a0 * PB1 - a1 * PB0 : 0.0;
            const double dBA0 = b0 > 0 ? -PA1 * (b0 - 1) + a1 * PB2 + at * PB1 : 0.0;
            const double dBA1 = b1 > 0 ? -PA2 * (b1 - 1) + at * PB0 + a0 * PB2 : 0.0;
            const double dBA2 = bt > 0 ? -PA0 * (bt - 1) + a0 * PB1 + a1 * PB0 : 0.0;
            push[2 * c0] = dAB0 - dAB2;
            push[2 * c0 + 1] = dBA0 - dBA2;
            push[2 * c1] = dAB1 - dAB2;
            push[2 * c1 + 1] = dBA1 - dBA2;
            G0 += dAB2;
            G1 += dBA2;
        } else {
            Agg0 h{};
            for (int k = 0; k < nc; ++k) agg0_add(h, av[kids[cb + k]], bv[kids[cb + k]], 1);
            agg0_add(h, at, bt, 1);
            const double f0 = agg0_f(h);
            fsum += f0;
            if (!want) continue;
            auto delta = [&](double a, double b, int dir) {
                Agg0 q = h;
                agg0_add(q, a, b, -1);
                if (dir == 0) agg0_add(q, a - 1, b + 1, 1); else agg0_add(q, a + 1, b - 1, 1);
                return agg0_f(q) - f0;
            };
            const double tAB = at > 0 ? delta(at, bt, 0) : 0.0, tBA = bt > 0 ? delta(at, bt, 1) : 0.0;
            for (int k = 0; k < nc; ++k) {
                const int c = kids[cb + k];
                push[2 * c] = (av[c] > 0 ? delta(av[c], bv[c], 0) : 0.0) - tAB;
                push[2 * c + 1] = (bv[c] > 0 ? delta(av[c], bv[c], 1) : 0.0) - tBA;
            }
            G0 += tAB;
            G1 += tBA;
        }
    }

    const double PA = (A * A - RA) / 2, PB = (B * B - RB) / 2, tq = PA * PB;
    const double score = 2 * fsum - tq;
    if (!want) return score;

    push[2 * root] = 0;
    push[2 * root + 1] = 0;
    const double dtot0 = (PA - (A - 1)) * (PB + B) - tq;
    const double dtot1 = (PA + A) * (PB - (B - 1)) - tq;
    const int* __restrict par = C_.parent.data();
    const int* __restrict leaf = C_.leaf_real.data();
    const int8_t* __restrict side = pr;
    for (int q = qe - 2; q >= qb; --q) {            // internal nodes, parents first (the root is last)
        const int v = inner[q], p = par[v];
        push[2 * v] += push[2 * p];
        push[2 * v + 1] += push[2 * p + 1];
    }
    const int* __restrict rl = rleaf_.data();
    for (int q = rleaf_off_[g]; q < rleaf_off_[g + 1]; ++q) {
        const int v = rl[q], p = par[v], i = leaf[v];
        real_gain[i] += side[i] == 0 ? 2 * ((push[2 * v] + push[2 * p]) + G0) - dtot0
                                     : 2 * ((push[2 * v + 1] + push[2 * p + 1]) + G1) - dtot1;
    }
    return score;
}

namespace {

// Binary node, branches 0, 1 = children and 2 = top: PA[i] pairs branches (i, i+1), f = sum_i PA[i] PB[i+2]
inline double binary_f(const double a[3], const double b[3], const double sA[3], const double sB[3],
                       double PA[3], double PB[3]) {
    for (int i = 0; i < 3; ++i) {
        PA[i] = a[i] * a[(i + 1) % 3] - sA[i];
        PB[i] = (b[i] * b[i] - sB[i]) / 2;
    }
    return PA[0] * PB[2] + PA[1] * PB[0] + PA[2] * PB[1];
}

// real-taxon deltas of a binary node: pushed to its children, the top branch's summed over the tree
inline void binary_push(double* push, int c0, int c1, const double a[3], const double b[3], const int ra3[3],
                        const int rb3[3], const double PA[3], const double PB[3], double& G0, double& G1) {
    double dAB[3], dBA[3];
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3, l = (i + 2) % 3;
        // the pair opposite branch i is PA[j]; pairs (i,j) and (l,i) lose/gain a_j and a_l
        dAB[i] = ra3[i] > 0 ? PA[j] * b[i] - a[j] * PB[l] - a[l] * PB[j] : 0.0;
        dBA[i] = rb3[i] > 0 ? -PA[j] * (b[i] - 1) + a[j] * PB[l] + a[l] * PB[j] : 0.0;
    }
    push[2 * c0] = dAB[0] - dAB[2];
    push[2 * c0 + 1] = dBA[0] - dBA[2];
    push[2 * c1] = dAB[1] - dAB[2];
    push[2 * c1 + 1] = dBA[1] - dBA[2];
    G0 += dAB[2];
    G1 += dBA[2];
}

// change of a binary node's f when a dummy taxon with branch weights xs changes side (s = -1: A -> B)
inline double binary_dummy_delta(const double a[3], const double b[3], const double sA[3], const double sB[3],
                                 const double xs[3], double s, double f) {
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
    return f2 - f;
}

}  // namespace

// Polytomy node v in a tree with dummy taxa. The dummies with weight below v are explicit, with local
// index q = their position in v's list; every other dummy of the tree lies in the top branch only and
// matters there just through the B-side squares. Returns f; with `want`, also pushes the real-taxon
// deltas, sets coef = (alpha, beta, gamma) for the absent dummies (see tree_pass_dummies) and leaves
// each listed dummy's delta in S.pdelta.
double Scorer::poly_node(int v, const double* X, double sqXB, double at, double bt, int rat, int rbt, bool want,
                         double& G0, double& G1, double coef[3], Scratch& S, const int8_t* pd) {
    const int cb = C_.child_off[v], nc = C_.child_off[v + 1] - cb;
    const int eb = xoff_[v], m = xoff_[v + 1] - eb;
    PolyNode& pn = S.pn;
    pn.begin(m);
    double sq_top = sqXB;                      // B-side squares of the dummies absent below v
    for (int q = 0; q < m; ++q) {
        const int j = xj_[eb + q];
        S.pos[j] = q;
        pn.inA[q] = pd[j] == 0;
        if (pd[j] != 0) sq_top -= X[j] * X[j];
    }
    for (int k = 0; k < nc; ++k) {
        const int c = C_.children[cb + k];
        pn.branch(a_[c], b_[c], rb_[c]);
        for (int t = xoff_[c]; t < xoff_[c + 1]; ++t) pn.entry(S.pos[xj_[t]], xw_[t]);
    }
    pn.branch(at, bt, rbt, sq_top);
    for (int q = 0; q < m; ++q) {
        const double x = X[xj_[eb + q]] - xw_[eb + q];
        if (x != 0) pn.entry(q, x);
    }
    const double f0 = pn.build();
    if (!want) return f0;

    S.gd.resize(2 * static_cast<size_t>(nc + 1));
    for (int k = 0; k < nc; ++k) {
        const int c = C_.children[cb + k];
        S.gd[2 * k] = ra_[c] > 0 ? pn.real_delta(k, 0) : 0.0;
        S.gd[2 * k + 1] = rb_[c] > 0 ? pn.real_delta(k, 1) : 0.0;
    }
    const double tAB = pn.real_delta(nc, 0), tBA = pn.real_delta(nc, 1);   // also needed without real taxa
    S.gd[2 * nc] = rat > 0 ? tAB : 0.0;
    S.gd[2 * nc + 1] = rbt > 0 ? tBA : 0.0;
    for (int k = 0; k < nc; ++k) {
        const int c = C_.children[cb + k];
        push_[2 * c] = S.gd[2 * k] - S.gd[2 * nc];
        push_[2 * c + 1] = S.gd[2 * k + 1] - S.gd[2 * nc + 1];
    }
    G0 += S.gd[2 * nc];
    G1 += S.gd[2 * nc + 1];
    const double gamma = pn.pairs_without(nc);
    coef[0] = tAB;
    coef[1] = tBA - gamma;
    coef[2] = gamma;
    S.pdelta.resize(m);
    for (int q = 0; q < m; ++q) S.pdelta[q] = pn.dummy_delta(q);
    return f0;
}

// tree_pass for a tree with dummy-taxon leaves. A node handles explicitly only the dummies with weight
// below it (its list). Any other dummy of the tree has all its weight X in the node's top branch, and
// switching its side changes f there by
//   alpha X (A -> B)  or  beta X + gamma X^2 (B -> A),
// where alpha and beta + gamma are the top branch's deltas for one real taxon (X = 1 behaves alike) and
// gamma sums PA over the pairs of the node's other branches; the dummy's own pairs vanish. The
// coefficients are summed over the tree's nodes once, and at the nodes above a dummy's leaves the
// polynomial is replaced by the exact change. Per tree: O(nodes + sum of list lengths), not O(nodes * d).
double Scorer::tree_pass_dummies(int g, bool want, double* real_gain, double* dummy_gain, Scratch& S,
                                 const int8_t* pr, const int8_t* pd) {
    const std::vector<int>& pres = present_[g];
    const double* X = &tx_[static_cast<size_t>(g) * d_];
    double sqA = 0, sqXB = 0;                  // squares of the dummies' weights in the tree, per side
    for (int j : pres) (pd[j] == 0 ? sqA : sqXB) += X[j] * X[j];
    double* ds = S.dsum.data();                // per dummy: summed exact-minus-polynomial changes
    if (want) for (int j : pres) ds[j] = 0;

    const int root = C_.tree_off[g + 1] - 1;
    const double A = ta_[g], B = tb_[g];
    const int RA = tra_[g], RB = trb_[g];
    double fsum = 0, G0 = 0, G1 = 0, Salpha = 0, Sbeta = 0, Sgamma = 0;
    const int qb = inner_off_[g], qe = inner_off_[g + 1];
    for (int q = qb; q < qe; ++q) {
        const int v = inner_[q];
        const int cb = C_.child_off[v], nc = C_.child_off[v + 1] - cb;
        const int eb = xoff_[v], ee = xoff_[v + 1];
        const double at = A - a_[v], bt = B - b_[v];
        const int rat = RA - ra_[v], rbt = RB - rb_[v];
        if (nc == 2) {
            // branches 0, 1 = children, 2 = top; sA/sB exclude same-dummy pairs (absent dummies only add
            // X^2 to the top branch's B squares)
            const int c0 = C_.children[cb], c1 = C_.children[cb + 1];
            const double a[3] = {a_[c0], a_[c1], at}, b[3] = {b_[c0], b_[c1], bt};
            const int ra3[3] = {ra_[c0], ra_[c1], rat}, rb3[3] = {rb_[c0], rb_[c1], rbt};
            double sA[3] = {0, 0, 0}, sB[3] = {double(rb3[0]), double(rb3[1]), double(rb3[2]) + sqXB};
            for (int e = eb; e < ee; ++e) {
                const int j = xj_[e];
                const double x0 = xc0_[e], x1 = xc1_[e], x2 = X[j] - xw_[e];
                if (pd[j] == 0) { sA[0] += x0 * x1; sA[1] += x1 * x2; sA[2] += x2 * x0; }
                else { sB[0] += x0 * x0; sB[1] += x1 * x1; sB[2] += x2 * x2 - X[j] * X[j]; }
            }
            double PA[3], PB[3];
            const double f = binary_f(a, b, sA, sB, PA, PB);
            fsum += f;
            if (!want) continue;
            binary_push(push_.data(), c0, c1, a, b, ra3, rb3, PA, PB, G0, G1);
            const double alpha = PA[0] * b[2] - a[0] * PB[1] - a[1] * PB[0];   // beta = -alpha, gamma = PA[0]
            Salpha += alpha;
            Sbeta -= alpha;
            Sgamma += PA[0];
            for (int e = eb; e < ee; ++e) {
                const int j = xj_[e];
                const bool in_a = pd[j] == 0;
                const double x = X[j], xs[3] = {xc0_[e], xc1_[e], x - xw_[e]};
                const double poly = in_a ? alpha * x : -alpha * x + PA[0] * x * x;
                ds[j] += binary_dummy_delta(a, b, sA, sB, xs, in_a ? -1.0 : 1.0, f) - poly;
            }
        } else {
            double coef[3];
            const double f0 = poly_node(v, X, sqXB, at, bt, rat, rbt, want, G0, G1, coef, S, pd);
            fsum += f0;
            if (!want) continue;
            Salpha += coef[0];
            Sbeta += coef[1];
            Sgamma += coef[2];
            for (int k = 0; k < ee - eb; ++k) {
                const int j = xj_[eb + k];
                const double x = X[j];
                const double poly = pd[j] == 0 ? coef[0] * x : coef[1] * x + coef[2] * x * x;
                ds[j] += S.pdelta[k] - poly;
            }
        }
    }

    // per-tree totals: w(PA(g)) w(PB(g)) = all quartets with two taxa on each side (Section 2.5.2)
    const double PA = (A * A - RA - sqA) / 2, PB = (B * B - RB - sqXB) / 2, tq = PA * PB;
    const double score = 2 * fsum - tq;
    if (!want) return score;

    // real taxa: gain = path sum of pushed deltas + top-branch deltas of every node
    push_[2 * root] = 0;
    push_[2 * root + 1] = 0;
    const double dtot0 = (PA - (A - 1)) * (PB + B) - tq;   // one real taxon A -> B
    const double dtot1 = (PA + A) * (PB - (B - 1)) - tq;   // B -> A
    for (int q = qe - 2; q >= qb; --q) {                   // internal nodes, parents first (the root is last)
        const int v = inner_[q], p = C_.parent[v];
        push_[2 * v] += push_[2 * p];
        push_[2 * v + 1] += push_[2 * p + 1];
    }
    for (int q = rleaf_off_[g]; q < rleaf_off_[g + 1]; ++q) {
        const int v = rleaf_[q], p = C_.parent[v], i = C_.leaf_real[v];
        real_gain[i] += pr[i] == 0 ? 2 * ((push_[2 * v] + push_[2 * p]) + G0) - dtot0
                                        : 2 * ((push_[2 * v + 1] + push_[2 * p + 1]) + G1) - dtot1;
    }
    for (int j : pres) {
        const double x = X[j];
        const bool in_a = pd[j] == 0;
        const double dsum = ds[j] + (in_a ? Salpha * x : Sbeta * x + Sgamma * x * x);
        const double dt = in_a ? (PA - x * (A - x)) * (PB + x * B) - tq : (PA + x * A) * (PB - x * (B - x)) - tq;
        dummy_gain[j] += 2 * dsum - dt;
    }
    return score;
}

// `want` is explicit: with no real (or no dummy) taxa the gain vectors are empty and data() may be null
double Scorer::tree_pass(int g, bool want, double* real_gain, double* dummy_gain, Scratch& S,
                         const int8_t* pr, const int8_t* pd) {
    return present_[g].empty() ? tree_pass_plain(g, want, real_gain, pr)
                               : tree_pass_dummies(g, want, real_gain, dummy_gain, S, pr, pd);
}

// A move only records itself: the next pass over the trees applies it to each chunk's trees just before
// scoring them, so the update is parallel and needs no barrier of its own. Moves are applied in order,
// so every total sees the same sequence of additions as when applied at once.
void Scorer::swap_real(int i) {
    moves_.push_back({i, false, part_r_[i]});
    part_r_[i] = static_cast<int8_t>(1 - part_r_[i]);
}

void Scorer::swap_dummy(int j) {
    moves_.push_back({j, true, part_d_[j]});
    part_d_[j] = static_cast<int8_t>(1 - part_d_[j]);
}

void Scorer::apply_moves(int c) {
    for (const Move& m : moves_) m.dummy ? move_dummy(c, m.unit, m.from) : move_real(c, m.unit, m.from);
}

void Scorer::move_real(int c, int i, int from) {
    for (int g = chunk_[c]; g < chunk_[c + 1]; ++g) {
        int v = C_.real_leaf[static_cast<size_t>(g) * R_ + i];
        if (v < 0) continue;
        for (; v >= 0; v = C_.parent[v]) {
            if (from == 0) { a_[v] -= 1; b_[v] += 1; --ra_[v]; ++rb_[v]; }
            else { a_[v] += 1; b_[v] -= 1; ++ra_[v]; --rb_[v]; }
        }
        if (from == 0) { ta_[g] -= 1; tb_[g] += 1; --tra_[g]; ++trb_[g]; }
        else { ta_[g] += 1; tb_[g] -= 1; ++tra_[g]; --trb_[g]; }
    }
}

void Scorer::move_dummy(int c, int j, int from) {
    const int g0 = chunk_[c], g1 = chunk_[c + 1];
    const int* nb = dnode_.data() + doff_[j];
    const int* ne = dnode_.data() + doff_[j + 1];
    const int vb = C_.tree_off[g0], ve = C_.tree_off[g1];
    for (const int* it = std::lower_bound(nb, ne, vb); it < ne && *it < ve; ++it) {
        const int v = *it;
        const double x = dw_[it - dnode_.data()];
        if (from == 0) { a_[v] -= x; b_[v] += x; } else { a_[v] += x; b_[v] -= x; }
    }
    for (int g = g0; g < g1; ++g) {
        const double x = tx_[static_cast<size_t>(g) * d_ + j];
        if (from == 0) { ta_[g] -= x; tb_[g] += x; } else { ta_[g] += x; tb_[g] -= x; }
    }
}

}  // namespace wqfm
