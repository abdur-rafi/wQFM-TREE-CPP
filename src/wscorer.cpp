#include "wscorer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wqfm {

WScorer::WScorer(const Forest& f, bool compress, const Options& o) : fb_(f, compress), opt_(o) {
    fb_.set_forest_only(true);
    fb_.set_length_weighting(o.length);
}

void WScorer::set_kappa(double kappa) {
    if (kappa != 1.0) throw std::runtime_error("weighted quartets: only kappa = 1 is supported");
}

void WScorer::set_dummy_gamma(double gamma) {
    if (gamma != 0.0) throw std::runtime_error("weighted quartets: only dummy gamma = 0 is supported");
}

void WScorer::setup(const Subproblem& sp, const DummyPool& pool) {
    fb_.setup(sp, pool);
    build();
}

void WScorer::setup_child(const Subproblem& child, const IScorer& parent, int side) {
    fb_.setup_child(child, parent, side);
    build();
}

void WScorer::set_pool(ThreadPool* pool) {
    pool_ = pool;
    ensure_scratch();
}

void WScorer::ensure_scratch() {
    const size_t n = pool_ ? static_cast<size_t>(pool_->threads()) : 1;
    if (scratch_.size() < n) scratch_.resize(n);
    const size_t N = C_ ? static_cast<size_t>(C_->nodes()) : 0;
    for (Scratch& sc : scratch_) if (sc.ovr.size() != N) sc.ovr.assign(N, -1);
}

void WScorer::run_chunks(const std::function<void(int, int)>& fn) {
    const int C = static_cast<int>(chunk_.size()) - 1;
    if (pool_ && C > 1) pool_->parallel_for(C, fn);
    else for (int c = 0; c < C; ++c) fn(c, 0);
}

// Per-subproblem tables from the forest just built.
void WScorer::build() {
    C_ = &fb_.forest();
    const SubForest& C = *C_;
    R_ = fb_.reals();
    d_ = fb_.dummies();
    const int N = C.nodes(), K = C.k;
    part_r_.assign(R_, 0);
    part_d_.assign(d_, 0);

    // branch factors; the two branches of a degree-2 root are one branch: its support goes to the
    // first child's branch, the second child's branch carries none (the lengths add up anyway)
    ef_.assign(N, 1.0);
    df_.assign(N, 0.0);
    auto sfactor = [&](double s) {                 // 1 - rescaled support; a branch without support: 0
        if (!opt_.support || std::isnan(s)) return 0.0;
        const double t = std::min(1.0, std::max(0.0, (s - opt_.smin) / (opt_.smax - opt_.smin)));
        return 1.0 - t;
    };
    for (int v = 0; v < N; ++v) {
        if (C.parent[v] < 0) continue;
        if (opt_.length && !std::isnan(C.blen[v])) ef_[v] = std::exp(-C.blen[v]);
        df_[v] = sfactor(C.supp[v]);
    }
    cpos_.assign(N, 0);
    for (int v = 0; v < N; ++v)
        for (int k = C.child_off[v]; k < C.child_off[v + 1]; ++k) cpos_[C.children[k]] = k - C.child_off[v];
    for (int g = 0; g < K; ++g) {
        const int root = C.tree_off[g + 1] - 1;
        if (C.child_off[root + 1] - C.child_off[root] != 2) continue;
        const int c0 = C.children[C.child_off[root]], c1 = C.children[C.child_off[root] + 1];
        df_[c0] = sfactor(combine_support(C.supp[c0], C.supp[c1]));
        df_[c1] = opt_.support ? 1.0 : 0.0;
    }

    // dummies below each node, with their weight discounted to the node (children before parents)
    dl_off_.assign(N + 1, 0);
    dl_k_.clear();
    dl_x_.clear();
    {
        std::vector<double> acc(d_, 0.0);
        std::vector<int> seen;
        for (int v = 0; v < N; ++v) {
            if (C.leaf_dummy[v] >= 0) {
                dl_k_.push_back(C.leaf_dummy[v]);
                dl_x_.push_back(C.leaf_w[v]);
            } else {
                seen.clear();
                for (int e = C.child_off[v]; e < C.child_off[v + 1]; ++e) {
                    const int c = C.children[e];
                    for (int t = dl_off_[c]; t < dl_off_[c + 1]; ++t) {
                        if (acc[dl_k_[t]] == 0) seen.push_back(dl_k_[t]);
                        acc[dl_k_[t]] += ef_[c] * dl_x_[t];
                    }
                }
                std::sort(seen.begin(), seen.end());
                for (int j : seen) { dl_k_.push_back(j); dl_x_.push_back(acc[j]); acc[j] = 0; }
            }
            dl_off_[v + 1] = static_cast<int>(dl_k_.size());
        }
    }
    // per tree its dummies (the root's list); per list entry the dummy's position among them; per
    // internal node the discounted weight of each of the tree's dummies in its top branch:
    //   xtop(v) = e_v (xtop(p) + x(p) - e_v x(v)) for parent p, 0 at the root
    pres_off_.assign(1, 0);
    pres_.clear();
    dl_t_.assign(dl_k_.size(), -1);
    xtop_off_.assign(N + 1, 0);
    xtop_.clear();
    std::vector<int> tpos(d_, -1);
    for (int g = 0; g < K; ++g) {
        const int b = C.tree_off[g], e = C.tree_off[g + 1], root = e - 1;
        for (int t = dl_off_[root]; t < dl_off_[root + 1]; ++t) { tpos[dl_k_[t]] = static_cast<int>(pres_.size() - pres_off_[g]); pres_.push_back(dl_k_[t]); }
        const int dt = static_cast<int>(pres_.size()) - pres_off_[g];
        pres_off_.push_back(static_cast<int>(pres_.size()));
        for (int v = b; v < e; ++v) {
            for (int t = dl_off_[v]; t < dl_off_[v + 1]; ++t) dl_t_[t] = tpos[dl_k_[t]];
            xtop_off_[v + 1] = xtop_off_[v] + (C.child_off[v + 1] > C.child_off[v] ? dt : 0);
        }
        xtop_.resize(xtop_off_[e], 0.0);
        std::vector<double> xv(dt), xp(dt);
        for (int v = e - 1; v >= b; --v) {
            if (C.child_off[v + 1] == C.child_off[v]) continue;
            const int p = C.parent[v];
            double* out = &xtop_[xtop_off_[v]];
            if (p < 0) { std::fill(out, out + dt, 0.0); continue; }
            std::fill(xv.begin(), xv.end(), 0.0);
            std::fill(xp.begin(), xp.end(), 0.0);
            for (int t = dl_off_[v]; t < dl_off_[v + 1]; ++t) xv[dl_t_[t]] = dl_x_[t];
            for (int t = dl_off_[p]; t < dl_off_[p + 1]; ++t) xp[dl_t_[t]] = dl_x_[t];
            const double* top = &xtop_[xtop_off_[p]];
            const double ev = ef_[v];
            for (int q = 0; q < dt; ++q) out[q] = ev * (top[q] + xp[q] - ev * xv[q]);
        }
        for (int t = dl_off_[root]; t < dl_off_[root + 1]; ++t) tpos[dl_k_[t]] = -1;
    }

    // per (tree, dummy) entry of pres_: the dummy's cluster leaves in that tree
    cl_off_.assign(pres_.size() + 1, 0);
    cl_.clear();
    for (int g = 0; g < K; ++g) {
        const int pb = pres_off_[g], dt = pres_off_[g + 1] - pb;
        for (int q = 0; q < dt; ++q) tpos[pres_[pb + q]] = q;
        std::vector<std::vector<int>> lists(dt);
        for (int v = C.tree_off[g]; v < C.tree_off[g + 1]; ++v)
            if (C.leaf_dummy[v] >= 0) lists[tpos[C.leaf_dummy[v]]].push_back(v);
        for (int q = 0; q < dt; ++q) {
            cl_.insert(cl_.end(), lists[q].begin(), lists[q].end());
            cl_off_[pb + q + 1] = static_cast<int>(cl_.size());
            tpos[pres_[pb + q]] = -1;
        }
    }
    nv_.assign(N, NodeVals{});
    mk_.assign(dl_k_.size(), 0.0); kmk_.assign(dl_k_.size(), 0.0);
    pk_off_.assign(N + 1, 0);
    for (int v = 0; v < N; ++v) {                     // pairs of an A and a B dummy: at most L*L/4
        const size_t L = dl_off_[v + 1] - dl_off_[v];
        pk_off_[v + 1] = pk_off_[v] + L * L / 4;
    }
    mkk_.assign(pk_off_[N], 0.0); kmkk_.assign(pk_off_[N], 0.0);
    rAtop_.assign(N, 0.0); rBtop_.assign(N, 0.0); aTop_.assign(N, 0.0); bTop_.assign(N, 0.0);
    c_.assign(N, 0.0); pp_.assign(N, 0.0); ga_.assign(N, 0.0); gb_.assign(N, 0.0); GA_.assign(N, 0.0); GB_.assign(N, 0.0);
    for (Scratch& sc : scratch_) sc.ovr.assign(N, -1);
    ge_.assign(static_cast<size_t>(N) * 8, 0.0);
    gmk_off_.assign(N + 1, 0);
    for (int v = 0; v < N; ++v)
        gmk_off_[v + 1] = gmk_off_[v] + (dl_off_[v + 1] - dl_off_[v]) * (C.child_off[v + 1] - C.child_off[v]);
    gmk_.assign(gmk_off_[N], 0.0); gkmk_.assign(gmk_off_[N], 0.0);
    for (auto* r : {&rPA_, &rKA_, &rPB_, &rKB_, &rM_, &rKM_, &H_}) r->assign(N, 0.0);
    rMk_.assign(dl_k_.size(), 0.0); rKMk_.assign(dl_k_.size(), 0.0);

    scale_ = 0;
    for (int g = 0; g < K; ++g) {
        double w = 0;
        for (int v = C.tree_off[g]; v < C.tree_off[g + 1]; ++v)
            if (C.leaf_real[v] >= 0 || C.leaf_dummy[v] >= 0) w += C.leaf_w[v];
        scale_ += w * w * w * w;
    }
    // chunks of whole trees, balanced by node count, fixed by the forest
    constexpr int kChunkNodes = 128, kMaxChunks = 64;
    const int Cn = std::max(1, std::min(kMaxChunks, N / kChunkNodes));
    chunk_.assign(1, 0);
    for (int c = 1; c < Cn; ++c) {
        const int target = static_cast<int>(static_cast<long long>(N) * c / Cn);
        int g = static_cast<int>(std::lower_bound(C.tree_off.begin(), C.tree_off.end(), target) - C.tree_off.begin());
        g = std::max(g, chunk_.back());
        if (g > chunk_.back() && g < K) chunk_.push_back(g);
    }
    chunk_.push_back(K);
    chunk_score_.assign(chunk_.size() - 1, 0.0);
    chunk_gain_.assign((chunk_.size() - 1) * static_cast<size_t>(R_ + d_), 0.0);
    ensure_scratch();
}

namespace {
inline void ensure(std::vector<double>& v, size_t n) { if (v.size() < n) v.resize(n); }
// the per-pair accumulators of Scratch::pacc, over the child branches i holding both dummies p and q:
// sum of x_p x_q; of M = Mkk(p, q) and M^2, K = KMkk(p, q) and K^2 (discounted to the node); and with
// w = M - K, of w, w x_p, w x_q, w x_p x_q
enum : int { kXX, kM, kM2, kK, kK2, kW, kWxp, kWxq, kWxx, kPair };
}  // namespace

// The branch values of node m: for each child branch i (then the top branch, index D) the values of
// the child, discounted to m: single sums by exp(-length), pair sums by 1 - support (K*) or not (P*, M).
// The child's dummy list is a sublist of m's. The scratch vectors only grow; entries are written here.
void WScorer::branch_values(int m, const int8_t* pd, Scratch& S) const {
    const SubForest& C = *C_;
    const int cb = C.child_off[m], D = C.child_off[m + 1] - cb;
    const int lb = dl_off_[m], L = dl_off_[m + 1] - lb;
    const int B = D + 1;
    S.D = D; S.L = L;
    for (auto* v : {&S.rA, &S.rB, &S.a, &S.b, &S.PA, &S.KA, &S.PB, &S.KB, &S.M, &S.KM}) ensure(*v, B);
    double *rA = S.rA.data(), *rB = S.rB.data(), *a = S.a.data(), *b = S.b.data();
    double *PA = S.PA.data(), *KA = S.KA.data(), *PB = S.PB.data(), *KB = S.KB.data(), *M = S.M.data(), *KM = S.KM.data();
    double *x = nullptr, *Mk = nullptr, *KMk = nullptr, *pacc = nullptr;
    if (L > 0) {
        const size_t nLB = static_cast<size_t>(L) * B, nLD = static_cast<size_t>(L) * D;
        ensure(S.x, nLB); ensure(S.Mk, nLD); ensure(S.KMk, nLD);
        x = S.x.data(); Mk = S.Mk.data(); KMk = S.KMk.data();
        std::fill_n(x, nLB, 0.0); std::fill_n(Mk, nLD, 0.0); std::fill_n(KMk, nLD, 0.0);
        if (S.pos.size() < static_cast<size_t>(L)) { S.pos.resize(L); S.rk.resize(L); S.cA.resize(L); S.cB.resize(L); }
        int na = 0, nb = 0;                                   // ranks of the node's dummies among their side
        for (int q = 0; q < L; ++q) S.rk[q] = pd[dl_k_[lb + q]] == 0 ? na++ : nb++;
        S.La = na; S.Lb = nb;
        const size_t nP = static_cast<size_t>(na) * nb * kPair;
        ensure(S.pacc, nP);
        pacc = S.pacc.data();
        std::fill_n(pacc, nP, 0.0);
    }
    for (int i = 0; i < D; ++i) {
        const int c = C.children[cb + i];
        const NodeVals& v = nv_[c];
        const double ec = ef_[c], dc = df_[c];
        rA[i] = ec * v.rA; rB[i] = ec * v.rB; a[i] = ec * v.a; b[i] = ec * v.b;
        PA[i] = v.PA; KA[i] = dc * v.KA; PB[i] = v.PB; KB[i] = dc * v.KB; M[i] = v.M; KM[i] = dc * v.KM;
        if (L == 0) continue;
        const int ob = dl_off_[c], Lc = dl_off_[c + 1] - ob;
        if (Lc == 0) continue;
        int* pos = S.pos.data();
        for (int t = 0, q = 0; t < Lc; ++t) {
            while (dl_k_[lb + q] != dl_k_[ob + t]) ++q;
            pos[t] = q;
            x[static_cast<size_t>(q) * B + i] = ec * dl_x_[ob + t];
            Mk[static_cast<size_t>(q) * D + i] = mk_[ob + t];
            KMk[static_cast<size_t>(q) * D + i] = dc * kmk_[ob + t];
        }
        // the child's pairs of an A and a B dummy (its stored block: A rows, B columns, in list order)
        int nA = 0, nB = 0;
        for (int t = 0; t < Lc; ++t) (pd[dl_k_[ob + t]] == 0 ? S.cA[nA++] : S.cB[nB++]) = t;
        for (int ia = 0; ia < nA; ++ia) {
            const int t1 = S.cA[ia];
            const double xp = ec * dl_x_[ob + t1];
            const double* vm = mkk_.data() + pk_off_[c] + static_cast<size_t>(ia) * nB;
            const double* vk = kmkk_.data() + pk_off_[c] + static_cast<size_t>(ia) * nB;
            double* row = pacc + static_cast<size_t>(S.rk[pos[t1]]) * S.Lb * kPair;
            for (int ib = 0; ib < nB; ++ib) {
                const int t2 = S.cB[ib];
                const double xq = ec * dl_x_[ob + t2], m = vm[ib], k = dc * vk[ib], w = m - k;
                double* acc = row + static_cast<size_t>(S.rk[pos[t2]]) * kPair;
                acc[kXX] += xp * xq;
                acc[kM] += m; acc[kM2] += m * m; acc[kK] += k; acc[kK2] += k * k;
                acc[kW] += w; acc[kWxp] += w * xp; acc[kWxq] += w * xq; acc[kWxx] += w * xp * xq;
            }
        }
    }
    rA[D] = rB[D] = 0;
    a[D] = S.useTop ? S.aT : aTop_[m];
    b[D] = S.useTop ? S.bT : bTop_[m];
    PA[D] = KA[D] = PB[D] = KB[D] = M[D] = KM[D] = 0;
    if (L > 0) {
        const double* top = &xtop_[xtop_off_[m]];
        for (int q = 0; q < L; ++q) x[static_cast<size_t>(q) * B + D] = top[dl_t_[lb + q]];
    }
}

// The bottom-up values of node m from its branch values (children only): own pairs anchored at m
// (members in different child branches, excluding pairs under one dummy) plus the children's.
void WScorer::node_values(int m, const int8_t* pd, Scratch& S, NodeVals& out) {
    const int D = S.D, L = S.L, B = D + 1, lb = dl_off_[m];
    const double *a = S.a.data(), *b = S.b.data();
    double A = 0, Bs = 0, Saa = 0, Sbb = 0, Sab = 0, rA = 0, rB = 0;
    double SPA = 0, SKA = 0, SPB = 0, SKB = 0, SM = 0, SKM = 0;
    for (int i = 0; i < D; ++i) {
        A += a[i]; Bs += b[i]; Saa += a[i] * a[i]; Sbb += b[i] * b[i]; Sab += a[i] * b[i];
        rA += S.rA[i]; rB += S.rB[i];
        SPA += S.PA[i]; SKA += S.KA[i]; SPB += S.PB[i]; SKB += S.KB[i]; SM += S.M[i]; SKM += S.KM[i];
    }
    double ownPA = (A * A - Saa) / 2, ownPB = (Bs * Bs - Sbb) / 2;
    const double ownM = A * Bs - Sab;
    out.rA = rA; out.rB = rB; out.a = rA; out.b = rB;
    if (L > 0) {
        const double *x = S.x.data(), *Mk = S.Mk.data(), *KMk = S.KMk.data();
        ensure(S.XC, L);
        double* XC = S.XC.data();
        for (int p = 0; p < L; ++p) {
            const double* xp = x + static_cast<size_t>(p) * B;
            double X = 0, sx2 = 0, sxb = 0, sax = 0;
            for (int i = 0; i < D; ++i) { X += xp[i]; sx2 += xp[i] * xp[i]; sxb += xp[i] * b[i]; sax += a[i] * xp[i]; }
            XC[p] = X;
            double own;
            if (pd[dl_k_[lb + p]] == 0) { ownPA -= (X * X - sx2) / 2; own = X * Bs - sxb; out.a += X; }
            else { ownPB -= (X * X - sx2) / 2; own = A * X - sax; out.b += X; }
            double sm = 0, skm = 0;
            for (int i = 0; i < D; ++i) { sm += Mk[static_cast<size_t>(p) * D + i]; skm += KMk[static_cast<size_t>(p) * D + i]; }
            mk_[lb + p] = own + sm;
            kmk_[lb + p] = own + skm;
        }
        // pairs (p in A, q in B): own ones split between child branches, plus the children's; stored as
        // the A x B block (see NodeVals), every entry written
        const double* pacc = S.pacc.data();
        const int Lb = S.Lb;
        for (int p = 0; p < L; ++p) {
            if (pd[dl_k_[lb + p]] != 0) continue;
            const size_t row = static_cast<size_t>(S.rk[p]) * Lb;
            for (int q = 0; q < L; ++q) {
                if (pd[dl_k_[lb + q]] != 1) continue;
                const size_t at = row + S.rk[q];
                const double* acc = pacc + at * kPair;
                const double own = XC[p] * XC[q] - acc[kXX];
                mkk_[pk_off_[m] + at] = own + acc[kM];
                kmkk_[pk_off_[m] + at] = own + acc[kK];
            }
        }
    }
    out.PA = ownPA + SPA; out.KA = ownPA + SKA;
    out.PB = ownPB + SPB; out.KB = ownPB + SKB;
    out.M = ownM + SM; out.KM = ownM + SKM;
}

// The quartets counted at m, from its branch values (children and top). With `grad`, also the
// gradients of the result with respect to each branch's values (it is multilinear: linear in the
// values of each branch), per branch i: a, b (every branch, the top included), and PA KA PB KB M KM
// and Mk KMk per dummy (child branches).
double WScorer::node_contrib(int m, const int8_t* pd, Scratch& S, bool grad, bool pairs) const {
    const int D = S.D, L = S.L, B = D + 1, lb = dl_off_[m];
    const double *a = S.a.data(), *b = S.b.data(), *PA = S.PA.data(), *KA = S.KA.data(), *PB = S.PB.data(), *KB = S.KB.data();
    const double *M = S.M.data(), *KM = S.KM.data(), *x = S.x.data(), *Mk = S.Mk.data(), *KMk = S.KMk.data();
    auto inA = [&](int p) { return pd[dl_k_[lb + p]] == 0; };
    auto inB = [&](int p) { return pd[dl_k_[lb + p]] == 1; };
    double Aall = 0, Ball = 0, Saa = 0, Sbb = 0, Sab = 0;
    for (int i = 0; i < B; ++i) { Aall += a[i]; Ball += b[i]; Saa += a[i] * a[i]; Sbb += b[i] * b[i]; Sab += a[i] * b[i]; }
    // pairs over all branches (members in different branches, not under one dummy); per dummy its
    // weight over all branches and its cross sum with the other side's single sums
    double Ahat = (Aall * Aall - Saa) / 2, Bhat = (Ball * Ball - Sbb) / 2;
    ensure(S.Xall, L); ensure(S.Sxb, L); ensure(S.Sxy, static_cast<size_t>(L) * L);
    ensure(S.SMk, L); ensure(S.SKMk, L);
    double *Xall = S.Xall.data(), *Sxb = S.Sxb.data(), *Sxy = S.Sxy.data(), *SMk = S.SMk.data(), *SKMk = S.SKMk.data();
    for (int p = 0; p < L; ++p) {
        const double* xp = x + static_cast<size_t>(p) * B;
        double X = 0, sx2 = 0, sv = 0;
        if (inA(p)) for (int i = 0; i < B; ++i) { X += xp[i]; sx2 += xp[i] * xp[i]; sv += xp[i] * b[i]; }
        else for (int i = 0; i < B; ++i) { X += xp[i]; sx2 += xp[i] * xp[i]; sv += a[i] * xp[i]; }
        Xall[p] = X; Sxb[p] = sv;
        (inA(p) ? Ahat : Bhat) -= (X * X - sx2) / 2;
    }
    (void)Sxy;
    // per branch i: the A-pairs and B-pairs among the other branches (An, Bn), and the mixed pairs
    // split among the other branches (Anp)
    ensure(S.An, B); ensure(S.Bn, B); ensure(S.Anp, B);
    double *An = S.An.data(), *Bn = S.Bn.data(), *Anp = S.Anp.data();
    for (int i = 0; i < B; ++i) {
        double an = Ahat - a[i] * (Aall - a[i]), bn = Bhat - b[i] * (Ball - b[i]);
        for (int p = 0; p < L; ++p) {
            const double xv = x[static_cast<size_t>(p) * B + i];
            (inA(p) ? an : bn) += xv * (Xall[p] - xv);
        }
        An[i] = an; Bn[i] = bn;
        Anp[i] = (Aall - a[i]) * (Ball - b[i]) - (Sab - a[i] * b[i]);
    }

    // satisfied: pairs in two child branches, or a pair in a child branch and the other split at m
    double SPA = 0, SKA = 0, SPB = 0, SKB = 0, SPAPB = 0, SKAKB = 0;
    for (int i = 0; i < D; ++i) {
        SPA += PA[i]; SKA += KA[i]; SPB += PB[i]; SKB += KB[i]; SPAPB += PA[i] * PB[i]; SKAKB += KA[i] * KB[i];
    }
    double c = SPA * SPB - SPAPB - (SKA * SKB - SKAKB);
    for (int i = 0; i < D; ++i) c += (PB[i] - KB[i]) * An[i] + (PA[i] - KA[i]) * Bn[i];
    // violated: mixed pairs in two child branches (inclusion-exclusion over the two sides' dummies),
    // or a mixed pair in a child branch and the other split at m
    double SM = 0, SM2 = 0, SKM = 0, SKM2 = 0;
    for (int i = 0; i < D; ++i) { SM += M[i]; SM2 += M[i] * M[i]; SKM += KM[i]; SKM2 += KM[i] * KM[i]; }
    double viol = (SM * SM - SM2) / 2 - (SKM * SKM - SKM2) / 2;
    for (int p = 0; p < L; ++p) {
        double sm = 0, sm2 = 0, k = 0, k2 = 0;
        for (int i = 0; i < D; ++i) {
            const double mv = Mk[static_cast<size_t>(p) * D + i], kv = KMk[static_cast<size_t>(p) * D + i];
            sm += mv; sm2 += mv * mv; k += kv; k2 += kv * kv;
        }
        SMk[p] = sm; SKMk[p] = k;
        viol -= (sm * sm - sm2) / 2 - (k * k - k2) / 2;
    }
    // terms over pairs of dummies (p in A, q in B), as in pair_term but from the per-pair accumulators:
    // with w_i = M_i - K_i, sum_i w_i ((X_p - x_pi)(X_q - x_qi) - (sxy - x_pi x_qi)) expands into them
    double pp = 0;
    if (pairs) {
        const double* pacc = S.pacc.data();
        for (int p = 0; p < L; ++p) {
            if (!inA(p)) continue;
            const double Xp = Xall[p], xpt = x[static_cast<size_t>(p) * B + D];
            const size_t row = static_cast<size_t>(S.rk[p]) * S.Lb;
            for (int q = 0; q < L; ++q) {
                if (!inB(q)) continue;
                const double* acc = pacc + (row + S.rk[q]) * kPair;
                const double Xq = Xall[q], sxy = acc[kXX] + xpt * x[static_cast<size_t>(q) * B + D];
                const double t = acc[kW] * (Xp * Xq - sxy) - Xp * acc[kWxq] - Xq * acc[kWxp] + 2 * acc[kWxx];
                pp += (acc[kM] * acc[kM] - acc[kM2]) / 2 - (acc[kK] * acc[kK] - acc[kK2]) / 2 + t;
            }
        }
    }
    S.pairPart = pp;
    viol += pp;
    ensure(S.E, static_cast<size_t>(L) * D);
    double* E = S.E.data();
    for (int i = 0; i < D; ++i) {
        viol += (M[i] - KM[i]) * Anp[i];
        for (int p = 0; p < L; ++p) {
            const double xv = x[static_cast<size_t>(p) * B + i], Xn = Xall[p] - xv;
            const double wp = Mk[static_cast<size_t>(p) * D + i] - KMk[static_cast<size_t>(p) * D + i];
            const double e = inA(p) ? Xn * (Ball - b[i]) - (Sxb[p] - xv * b[i]) : (Aall - a[i]) * Xn - (Sxb[p] - a[i] * xv);
            E[static_cast<size_t>(p) * D + i] = e;
            viol -= wp * e;
        }
    }
    const double result = c - viol;
    if (!grad) return result;

    // gradients w.r.t. a and b of branch i: sums over the child branches j != i, each formed once over
    // all child branches and corrected for j == i (linear in D; i == D, the top, excludes nothing)
    ensure(S.gA, B); ensure(S.gB, B);
    double U = 0, Ua = 0, V = 0, Vb = 0, Wt = 0, Wa = 0, Wb = 0;
    for (int j = 0; j < D; ++j) {
        const double u = PB[j] - KB[j], v = PA[j] - KA[j], w = M[j] - KM[j];
        U += u; Ua += u * a[j]; V += v; Vb += v * b[j]; Wt += w; Wa += w * a[j]; Wb += w * b[j];
    }
    for (int i = 0; i < B; ++i) {
        double u = 0, v = 0, w = 0;
        if (i < D) { u = PB[i] - KB[i]; v = PA[i] - KA[i]; w = M[i] - KM[i]; }
        S.gA[i] = (Aall - a[i]) * (U - u) - (Ua - u * a[i]) - ((Ball - b[i]) * (Wt - w) - (Wb - w * b[i]));
        S.gB[i] = (Ball - b[i]) * (V - v) - (Vb - v * b[i]) - ((Aall - a[i]) * (Wt - w) - (Wa - w * a[i]));
    }
    // per dummy p: sum over child branches j != i of (Mk - KMk)(p, j) times p's weight outside j and i
    for (int p = 0; p < L; ++p) {
        const double* xp = x + static_cast<size_t>(p) * B;
        const double* mkp = Mk + static_cast<size_t>(p) * D;
        const double* kmkp = KMk + static_cast<size_t>(p) * D;
        double WP = 0, WPX = 0;
        for (int j = 0; j < D; ++j) { const double wp = mkp[j] - kmkp[j]; WP += wp; WPX += wp * xp[j]; }
        double* g = inA(p) ? S.gB.data() : S.gA.data();
        for (int i = 0; i < B; ++i) {
            const double wp = i < D ? mkp[i] - kmkp[i] : 0.0;
            g[i] += (Xall[p] - xp[i]) * (WP - wp) - (WPX - wp * xp[i]);
        }
    }
    for (auto* v : {&S.gPA, &S.gKA, &S.gPB, &S.gKB, &S.gM, &S.gKM}) ensure(*v, D);
    for (int i = 0; i < D; ++i) {
        S.gPA[i] = (SPB - PB[i]) + Bn[i];
        S.gKA[i] = -(SKB - KB[i]) - Bn[i];
        S.gPB[i] = (SPA - PA[i]) + An[i];
        S.gKB[i] = -(SKA - KA[i]) - An[i];
        S.gM[i] = -((SM - M[i]) + Anp[i]);
        S.gKM[i] = (SKM - KM[i]) + Anp[i];
    }
    if (L > 0) {
        ensure(S.gMk, static_cast<size_t>(L) * D); ensure(S.gKMk, static_cast<size_t>(L) * D);
        for (int p = 0; p < L; ++p)
            for (int i = 0; i < D; ++i) {
                const size_t at = static_cast<size_t>(p) * D + i;
                S.gMk[at] = (SMk[p] - Mk[at]) + E[at];
                S.gKMk[at] = -((SKMk[p] - KMk[at]) + E[at]);
            }
    }
    return result;
}

// One tree: the discounted real weights below each node, then the top-branch sums (top-down), then
// bottom-up the values and contribution of every node (kept per node, with gradients and their
// subtree sums, when `keep`). Returns the tree's score.
double WScorer::tree_pass(int g, const int8_t* pr, const int8_t* pd, bool keep, Scratch& S) {
    const SubForest& C = *C_;
    const int b = C.tree_off[g], e = C.tree_off[g + 1];
    for (int v = b; v < e; ++v) {
        NodeVals& nv = nv_[v];
        const int cb = C.child_off[v], ce = C.child_off[v + 1];
        if (cb == ce) {
            nv.rA = nv.rB = 0;
            if (C.leaf_real[v] >= 0) (pr[C.leaf_real[v]] == 0 ? nv.rA : nv.rB) = 1;
            continue;
        }
        double rA = 0, rB = 0;
        for (int k = cb; k < ce; ++k) { const int c = C.children[k]; rA += ef_[c] * nv_[c].rA; rB += ef_[c] * nv_[c].rB; }
        nv.rA = rA; nv.rB = rB;
    }
    // top branch: r_top(v) = e_v (r_top(p) + r(p) - e_v r(v)), then the dummies on each side
    const int pb = pres_off_[g], dt = pres_off_[g + 1] - pb;
    for (int v = e - 1; v >= b; --v) {
        if (C.child_off[v + 1] == C.child_off[v]) continue;
        const int p = C.parent[v];
        if (p < 0) { rAtop_[v] = rBtop_[v] = 0; }
        else {
            const double ev = ef_[v];
            rAtop_[v] = ev * (rAtop_[p] + nv_[p].rA - ev * nv_[v].rA);
            rBtop_[v] = ev * (rBtop_[p] + nv_[p].rB - ev * nv_[v].rB);
        }
        double at = rAtop_[v], bt = rBtop_[v];
        const double* top = &xtop_[xtop_off_[v]];
        for (int q = 0; q < dt; ++q) (pd[pres_[pb + q]] == 0 ? at : bt) += top[q];
        aTop_[v] = at;
        bTop_[v] = bt;
    }
    double s = 0;
    for (int v = b; v < e; ++v) {
        NodeVals& nv = nv_[v];
        const int cb = C.child_off[v], D = C.child_off[v + 1] - cb;
        if (D == 0) {
            nv.a = nv.rA; nv.b = nv.rB;
            if (C.leaf_dummy[v] >= 0) (pd[C.leaf_dummy[v]] == 0 ? nv.a : nv.b) = C.leaf_w[v];
            nv.PA = nv.KA = nv.PB = nv.KB = nv.M = nv.KM = 0;
            for (int t = dl_off_[v]; t < dl_off_[v + 1]; ++t) mk_[t] = kmk_[t] = 0;   // a leaf has no pairs
            if (keep) GA_[v] = GB_[v] = 0;
            continue;
        }
        if (D == 2 && dl_off_[v + 1] == dl_off_[v]) {   // binary, no dummy below: node_values and
            const int c0 = C.children[cb], c1 = C.children[cb + 1];   // node_contrib in closed form
            const NodeVals& v0 = nv_[c0];
            const NodeVals& v1 = nv_[c1];
            const double e0 = ef_[c0], e1 = ef_[c1], d0 = df_[c0], d1 = df_[c1];
            const double a0 = e0 * v0.a, a1 = e1 * v1.a, a2 = aTop_[v];
            const double b0 = e0 * v0.b, b1 = e1 * v1.b, b2 = bTop_[v];
            const double PA0 = v0.PA, PA1 = v1.PA, KA0 = d0 * v0.KA, KA1 = d1 * v1.KA;
            const double PB0 = v0.PB, PB1 = v1.PB, KB0 = d0 * v0.KB, KB1 = d1 * v1.KB;
            const double M0 = v0.M, M1 = v1.M, KM0 = d0 * v0.KM, KM1 = d1 * v1.KM;
            const double oPA = a0 * a1, oPB = b0 * b1, oM = a0 * b1 + a1 * b0;   // pairs split at v
            nv.a = nv.rA; nv.b = nv.rB;
            nv.PA = oPA + PA0 + PA1; nv.KA = oPA + KA0 + KA1;
            nv.PB = oPB + PB0 + PB1; nv.KB = oPB + KB0 + KB1;
            nv.M = oM + M0 + M1; nv.KM = oM + KM0 + KM1;
            // per branch: A-pairs and B-pairs between the other two branches, and mixed pairs there
            const double An0 = a1 * a2, An1 = a0 * a2, Bn0 = b1 * b2, Bn1 = b0 * b2;
            const double Anp0 = a1 * b2 + a2 * b1, Anp1 = a0 * b2 + a2 * b0;
            const double w0 = M0 - KM0, w1 = M1 - KM1;
            const double cv = PA0 * PB1 + PA1 * PB0 - (KA0 * KB1 + KA1 * KB0)
                            + (PB0 - KB0) * An0 + (PA0 - KA0) * Bn0 + (PB1 - KB1) * An1 + (PA1 - KA1) * Bn1
                            - (M0 * M1 - KM0 * KM1 + w0 * Anp0 + w1 * Anp1);
            s += cv;
            if (!keep) continue;
            c_[v] = cv;
            pp_[v] = 0;
            double* g0 = &ge_[static_cast<size_t>(c0) * 8];
            double* g1 = &ge_[static_cast<size_t>(c1) * 8];
            g0[0] = (PB1 - KB1) * a2 - w1 * b2;
            g1[0] = (PB0 - KB0) * a2 - w0 * b2;
            g0[1] = (PA1 - KA1) * b2 - w1 * a2;
            g1[1] = (PA0 - KA0) * b2 - w0 * a2;
            g0[2] = PB1 + Bn0;  g1[2] = PB0 + Bn1;
            g0[3] = -KB1 - Bn0; g1[3] = -KB0 - Bn1;
            g0[4] = PA1 + An0;  g1[4] = PA0 + An1;
            g0[5] = -KA1 - An0; g1[5] = -KA0 - An1;
            g0[6] = -(M1 + Anp0); g1[6] = -(M0 + Anp1);
            g0[7] = KM1 + Anp0;   g1[7] = KM0 + Anp1;
            ga_[v] = (PB0 - KB0) * a1 - w0 * b1 + (PB1 - KB1) * a0 - w1 * b0;
            gb_[v] = (PA0 - KA0) * b1 - w0 * a1 + (PA1 - KA1) * b0 - w1 * a0;
            GA_[v] = ga_[v] + e0 * GA_[c0] + e1 * GA_[c1];
            GB_[v] = gb_[v] + e0 * GB_[c0] + e1 * GB_[c1];
            continue;
        }
        branch_values(v, pd, S);
        node_values(v, pd, S, nv);
        const double cv = node_contrib(v, pd, S, keep);
        s += cv;
        if (!keep) continue;
        c_[v] = cv;
        pp_[v] = S.pairPart;
        const int L = S.L;
        ga_[v] = S.gA[D]; gb_[v] = S.gB[D];
        double GA = ga_[v], GB = gb_[v];
        for (int i = 0; i < D; ++i) {
            const int c = C.children[cb + i];
            double* gr = &ge_[static_cast<size_t>(c) * 8];
            gr[0] = S.gA[i]; gr[1] = S.gB[i]; gr[2] = S.gPA[i]; gr[3] = S.gKA[i]; gr[4] = S.gPB[i]; gr[5] = S.gKB[i]; gr[6] = S.gM[i]; gr[7] = S.gKM[i];
            GA += ef_[c] * GA_[c]; GB += ef_[c] * GB_[c];
        }
        for (int p = 0; p < L; ++p)
            for (int i = 0; i < D; ++i) {
                gmk_[gmk_off_[v] + static_cast<size_t>(p) * D + i] = S.gMk[static_cast<size_t>(p) * D + i];
                gkmk_[gmk_off_[v] + static_cast<size_t>(p) * D + i] = S.gKMk[static_cast<size_t>(p) * D + i];
            }
        GA_[v] = GA; GB_[v] = GB;
    }
    return s;
}

// All real-taxon gains of tree g, after tree_pass with `keep`. For a taxon t moving from A to B, a
// node m above t changes by the gradients of its contribution w.r.t. the branch towards t times the
// changes of that branch's values: a and b by -+ e(t->m); the pair sums by the pairs t forms with the
// units joining it at each ancestor w below m, weighted by e(t->w) (and by the support path from m to
// w for the K sums); so the gain is sum over w of e(t->w) times a value at w that sums the gradients of
// all nodes above w (R, top-down). Nodes off the path see t in their top branch (GA, GB). One top-down
// pass accumulates the path sums H for every leaf; from B to A the gain is -H.
void WScorer::gain_pass(int g, const int8_t* pr, const int8_t* pd, double* real_gain, Scratch& S) {
    const SubForest& C = *C_;
    const int b = C.tree_off[g], e = C.tree_off[g + 1];
    // R: gradient sums over the ancestors, per node (internal nodes only), parents first
    for (int w = e - 1; w >= b; --w) {
        if (C.child_off[w + 1] == C.child_off[w]) continue;
        const int p = C.parent[w], lb = dl_off_[w], L = dl_off_[w + 1] - lb;
        if (p < 0) {
            rPA_[w] = rKA_[w] = rPB_[w] = rKB_[w] = rM_[w] = rKM_[w] = 0;
            for (int t = 0; t < L; ++t) rMk_[lb + t] = rKMk_[lb + t] = 0;
            H_[w] = 0;
            continue;
        }
        const double* gr = &ge_[static_cast<size_t>(w) * 8];
        const double dw = df_[w];
        rPA_[w] = gr[2] + rPA_[p]; rKA_[w] = dw * (gr[3] + rKA_[p]);
        rPB_[w] = gr[4] + rPB_[p]; rKB_[w] = dw * (gr[5] + rKB_[p]);
        rM_[w] = gr[6] + rM_[p]; rKM_[w] = dw * (gr[7] + rKM_[p]);
        const int iw = cpos_[w];
        const int Dp = C.child_off[p + 1] - C.child_off[p], pb = dl_off_[p];
        for (int t = 0, q = 0; t < L; ++t) {
            while (dl_k_[pb + q] != dl_k_[lb + t]) ++q;
            const size_t at = gmk_off_[p] + static_cast<size_t>(q) * Dp + iw;
            rMk_[lb + t] = gmk_[at] + rMk_[pb + q];
            rKMk_[lb + t] = dw * (gkmk_[at] + rKMk_[pb + q]);
        }
    }
    // H: for each child c of w, the value at w seen from c (the other branches' units joining there),
    // discounted down to c; H(c) = e_c (q(w, c) + H(w))
    for (int w = e - 1; w >= b; --w) {
        const int cb = C.child_off[w], D = C.child_off[w + 1] - cb;
        if (D == 0) continue;
        const int lb = dl_off_[w], L = dl_off_[w + 1] - lb;
        double tot = 0;
        for (int i = 0; i < D; ++i) { const int c = C.children[cb + i]; tot += ef_[c] * (GB_[c] - GA_[c]); }
        const double cA = -rPA_[w] - rKA_[w] + rM_[w] + rKM_[w], cB = rPB_[w] + rKB_[w] - rM_[w] - rKM_[w];
        for (int i = 0; i < D; ++i) {
            const int c = C.children[cb + i];
            const double ec = ef_[c];
            const double* gr = &ge_[static_cast<size_t>(c) * 8];
            double q = (nv_[w].a - ec * nv_[c].a) * cA + (nv_[w].b - ec * nv_[c].b) * cB
                     + (gr[1] - gr[0]) + (tot - ec * (GB_[c] - GA_[c]));
            const int ob = dl_off_[c], Lc = dl_off_[c + 1] - ob;
            for (int t = 0, u = 0; t < L; ++t) {
                double xc = 0;
                if (u < Lc && dl_k_[ob + u] == dl_k_[lb + t]) { xc = dl_x_[ob + u]; ++u; }
                const double beta = dl_x_[lb + t] - ec * xc;
                const double sgn = pd[dl_k_[lb + t]] == 0 ? 1.0 : -1.0;
                q += sgn * beta * (rMk_[lb + t] + rKMk_[lb + t]);
            }
            H_[c] = ec * (q + H_[w]);
        }
    }
    for (int v = b; v < e; ++v) {
        const int r = C.leaf_real[v];
        if (r >= 0) real_gain[r] += pr[r] == 0 ? H_[v] : -H_[v];
    }
    (void)S;
}

// The violated quartets with both A taxa under dummy p and both B taxa under dummy q that the
// inclusion-exclusion over the two sides counts twice, at the node of S's branch values: M, K are the
// pair sums of (p, q) per child branch (plain and support-discounted), and the geometry comes from x.
double WScorer::pair_term(int p, int q, const double* M, const double* K, Scratch& S) const {
    const int D = S.D, B = D + 1;
    const double* xp = &S.x[static_cast<size_t>(p) * B];
    const double* xq = &S.x[static_cast<size_t>(q) * B];
    double Xp = 0, Xq = 0, sxy = 0;
    for (int i = 0; i < B; ++i) { Xp += xp[i]; Xq += xq[i]; sxy += xp[i] * xq[i]; }
    double sm = 0, sm2 = 0, k = 0, k2 = 0, t = 0;
    for (int i = 0; i < D; ++i) {
        sm += M[i]; sm2 += M[i] * M[i]; k += K[i]; k2 += K[i] * K[i];
        t += (M[i] - K[i]) * ((Xp - xp[i]) * (Xq - xq[i]) - (sxy - xp[i] * xq[i]));
    }
    return (sm * sm - sm2) / 2 - (k * k - k2) / 2 + t;
}

// Change of tree g's score when dummy j (the tree's tq-th dummy) switches sides, after tree_pass with
// `keep`. Ancestors of j's clusters, children first: the branch values are loaded unchanged and
// adjusted for the switch from j's pair sums along the path, the contribution without the dummy-pair
// terms is recomputed, and of those terms only the ones with j are replaced. A subtree hanging off
// the ancestors sees j only in its top branch: it changes by j's discounted weight at its root times
// its gradient sum.
double WScorer::dummy_change(int g, int tq, int j, Scratch& S) {
    const SubForest& C = *C_;
    const int from = part_d_[j];
    const double sgn = from == 0 ? 1.0 : -1.0;
    const int qq = pres_off_[g] + tq;
    S.path.clear();
    S.pv.clear();
    S.W.clear();
    S.KW.clear();
    for (int e = cl_off_[qq]; e < cl_off_[qq + 1]; ++e) {
        const int leaf = cl_[e];
        S.ovr[leaf] = -3;                        // a cluster of j: no pairs below it
        for (int v = C.parent[leaf]; v >= 0 && S.ovr[v] == -1; v = C.parent[v]) { S.ovr[v] = -2; S.path.push_back(v); }
    }
    std::sort(S.path.begin(), S.path.end());
    double delta = 0;
    for (int p : S.path) {
        const int lb = dl_off_[p], L = dl_off_[p + 1] - lb, cb = C.child_off[p], D = C.child_off[p + 1] - cb;
        const int tk = static_cast<int>(std::lower_bound(dl_k_.begin() + lb, dl_k_.begin() + lb + L, j) - dl_k_.begin()) - lb;
        const double xk = dl_x_[lb + tk];                 // j's weight below p
        const double xt = xtop_[xtop_off_[p] + tq];       // j's weight in p's top branch
        // p's contribution after the switch, and j's pair sums at p for its parent (me)
        Scratch::PathVals me{0, 0, 0};
        delta += D == 2 ? dummy_change_binary(p, j, tk, from, xt, S, me) : dummy_change_poly(p, j, tk, from, xt, S, me);
        S.ovr[p] = static_cast<int>(S.pv.size());
        S.pv.push_back(me);
        // subtrees hanging off the path see j only in their top branch
        for (int k = cb; k < cb + D; ++k) {
            const int c = C.children[k];
            if (S.ovr[c] != -1 || C.child_off[c + 1] == C.child_off[c]) continue;
            delta += sgn * ef_[c] * (xt + xk) * (GB_[c] - GA_[c]);
        }
    }
    for (int e = cl_off_[qq]; e < cl_off_[qq + 1]; ++e) S.ovr[cl_[e]] = -1;
    for (int p : S.path) S.ovr[p] = -1;
    return delta;
}

// dummy_change at a binary node p on j's ancestor path: the change of p's contribution when j switches
// sides, and j's pair sums at p (me) for p's parent. Same arithmetic as the generic path (branch values
// adjusted for the switch, contribution without the dummy-pair terms, j's dummy-pair terms replaced),
// specialised to two children and fused into one pass over p's dummy list, which is the union of the
// children's lists (merged on the fly, so nothing is stored per dummy).
double WScorer::dummy_change_binary(int p, int j, int tk, int from, double xt, Scratch& S, Scratch::PathVals& me) {
    const SubForest& C = *C_;
    const int lb = dl_off_[p], L = dl_off_[p + 1] - lb, cb = C.child_off[p];
    const double sgn = from == 0 ? 1.0 : -1.0;
    const double* top = &xtop_[xtop_off_[p]];
    // the children: branch values, their dummy lists, j's position and weight there, and j's path values
    double a[3], b[3], PA[2], KA[2], PB[2], KB[2], M[2], KM[2], SS[2] = {0, 0}, KSS[2] = {0, 0}, OP[2] = {0, 0}, KOP[2] = {0, 0};
    double ec[2], dc[2], xj[2] = {0, 0};
    int st[2], ob[2], Lc[2], tjc[2] = {-1, -1}, rj[2] = {0, 0}, nbc[2] = {0, 0};
    size_t pkc[2];
    const Scratch::PathVals* cv[2] = {nullptr, nullptr};
    for (int i = 0; i < 2; ++i) {
        const int c = C.children[cb + i];
        const NodeVals& v = nv_[c];
        pkc[i] = pk_off_[c];
        ec[i] = ef_[c]; dc[i] = df_[c];
        st[i] = S.ovr[c];                           // -1: j not below c; -3: a cluster of j; >= 0: on j's path
        a[i] = ec[i] * v.a; b[i] = ec[i] * v.b;
        PA[i] = v.PA; KA[i] = dc[i] * v.KA; PB[i] = v.PB; KB[i] = dc[i] * v.KB; M[i] = v.M; KM[i] = dc[i] * v.KM;
        ob[i] = dl_off_[c]; Lc[i] = dl_off_[c + 1] - ob[i];
        if (st[i] != -1) {
            tjc[i] = static_cast<int>(std::lower_bound(dl_k_.begin() + ob[i], dl_k_.begin() + ob[i] + Lc[i], j) - dl_k_.begin()) - ob[i];
            xj[i] = ec[i] * dl_x_[ob[i] + tjc[i]];
        }
        if (st[i] >= 0) {
            cv[i] = &S.pv[st[i]];
            SS[i] = cv[i]->SS; KSS[i] = dc[i] * cv[i]->KSS; OP[i] = mk_[ob[i] + tjc[i]]; KOP[i] = dc[i] * kmk_[ob[i] + tjc[i]];
            int na = 0, nb = 0;                     // j's rank on its side, and the child's B dummies
            for (int t = 0; t < Lc[i]; ++t) {
                const bool inA = part_d_[dl_k_[ob[i] + t]] == 0;
                if (t == tjc[i]) rj[i] = inA ? na : nb;
                if (inA) ++na; else ++nb;
            }
            nbc[i] = nb;
        }
    }
    a[2] = aTop_[p] - sgn * xt;
    b[2] = bTop_[p] + sgn * xt;
    const double xj0 = xj[0], xj1 = xj[1], Xj = xj0 + xj1 + xt;
    {   // j's pairs with its own side, split at p (before the switch), for its parent
        const double s0 = from == 0 ? a[0] : b[0], s1 = from == 0 ? a[1] : b[1];
        const double own = xj0 * s1 + xj1 * s0 - 2 * xj0 * xj1;
        me.SS = own + SS[0] + SS[1];
        me.KSS = own + KSS[0] + KSS[1];
    }
    // the switch, in the children's single and pair sums (the per-dummy sums: in the pass below)
    for (int i = 0; i < 2; ++i) {
        if (st[i] == -1) continue;
        if (from == 0) { a[i] -= xj[i]; b[i] += xj[i]; } else { a[i] += xj[i]; b[i] -= xj[i]; }
        if (st[i] == -3) continue;
        if (from == 0) { PA[i] -= SS[i]; KA[i] -= KSS[i]; PB[i] += OP[i]; KB[i] += KOP[i]; }
        else { PB[i] -= SS[i]; KB[i] -= KSS[i]; PA[i] += OP[i]; KA[i] += KOP[i]; }
        M[i] += SS[i] - OP[i];
        KM[i] += KSS[i] - KOP[i];
    }
    const double Aall = a[0] + a[1] + a[2], Ball = b[0] + b[1] + b[2];
    const double Saa = a[0] * a[0] + a[1] * a[1] + a[2] * a[2], Sbb = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    const double Sab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    double Ahat = (Aall * Aall - Saa) / 2, Bhat = (Ball * Ball - Sbb) / 2;
    double exA0 = 0, exA1 = 0, exB0 = 0, exB1 = 0, mkboth = 0, eterm = 0, ppd = 0;
    // one pass over p's dummies: their weights per branch, the children's per-dummy sums (switched), j's
    // pair sums with each (for the parent), the contribution's per-dummy terms and the dummy-pair terms
    me.woff = static_cast<int>(S.W.size());
    S.W.resize(S.W.size() + L, 0.0);
    S.KW.resize(S.KW.size() + L, 0.0);
    double* W = &S.W[me.woff];
    double* KW = &S.KW[me.woff];
    int u[2] = {0, 0}, ra[2] = {0, 0}, rb[2] = {0, 0};
    for (int t = 0; t < L; ++t) {
        const int k = dl_k_[lb + t];
        const int before = part_d_[k];
        double y[2] = {0, 0}, m[2] = {0, 0}, km[2] = {0, 0}, w[2] = {0, 0}, kw[2] = {0, 0}, mo[2] = {0, 0}, kmo[2] = {0, 0};
        for (int i = 0; i < 2; ++i) {
            if (u[i] >= Lc[i] || dl_k_[ob[i] + u[i]] != k) continue;
            const int tc = u[i]++;
            const int r = before == 0 ? ra[i]++ : rb[i]++;  // its rank among the child's dummies of its side
            y[i] = ec[i] * dl_x_[ob[i] + tc];
            m[i] = mk_[ob[i] + tc];
            km[i] = dc[i] * kmk_[ob[i] + tc];
            if (!cv[i] || t == tk) continue;
            w[i] = S.W[cv[i]->woff + tc];
            kw[i] = dc[i] * S.KW[cv[i]->woff + tc];
            if (before != from) {                   // (j, t) mixed: j's stored pair sums with t
                const size_t sx = from == 0 ? static_cast<size_t>(rj[i]) * nbc[i] + r : static_cast<size_t>(r) * nbc[i] + rj[i];
                mo[i] = mkk_[pkc[i] + sx];
                kmo[i] = dc[i] * kmkk_[pkc[i] + sx];
            }
        }
        const double y0 = y[0], y1 = y[1], y2 = top[dl_t_[lb + t]], X = y0 + y1 + y2;
        const int side = t == tk ? 1 - from : before;
        if (t != tk) {
            const double own = y0 * xj1 + y1 * xj0;
            W[t] = own + w[0] + w[1];
            KW[t] = own + kw[0] + kw[1];
        }
        for (int i = 0; i < 2; ++i) {
            if (!cv[i]) continue;
            if (t == tk) { m[i] = SS[i]; km[i] = KSS[i]; }
            else {
                const double sg = before == from ? 1.0 : -1.0;     // same side as j before: (t, j) becomes mixed
                m[i] += sg * w[i];
                km[i] += sg * kw[i];
            }
        }
        // the contribution without the dummy-pair terms (node_contrib with two children), after the switch
        const double split = y0 * y1 + y0 * y2 + y1 * y2;           // the dummy's own pairs across branches
        const double w0 = m[0] - km[0], w1 = m[1] - km[1];
        mkboth += m[0] * m[1] - km[0] * km[1];
        if (side == 0) {
            Ahat -= split;
            exA0 += y0 * (X - y0); exA1 += y1 * (X - y1);
            const double sxb = y0 * b[0] + y1 * b[1] + y2 * b[2];
            eterm += w0 * ((X - y0) * (Ball - b[0]) - (sxb - y0 * b[0])) + w1 * ((X - y1) * (Ball - b[1]) - (sxb - y1 * b[1]));
        } else {
            Bhat -= split;
            exB0 += y0 * (X - y0); exB1 += y1 * (X - y1);
            const double sax = a[0] * y0 + a[1] * y1 + a[2] * y2;
            eterm += w0 * ((Aall - a[0]) * (X - y0) - (sax - a[0] * y0)) + w1 * ((Aall - a[1]) * (X - y1) - (sax - a[1] * y1));
        }
        // dummy-pair terms: j's old pairs out, its new ones in
        if (t == tk) continue;
        const double sxy = y0 * xj0 + y1 * xj1 + y2 * xt;
        const double g0 = (X - y0) * (Xj - xj0) - (sxy - y0 * xj0), g1 = (X - y1) * (Xj - xj1) - (sxy - y1 * xj1);
        if (before != from)                         // (j, t) was a mixed pair: remove its term
            ppd -= mo[0] * mo[1] - kmo[0] * kmo[1] + (mo[0] - kmo[0]) * g0 + (mo[1] - kmo[1]) * g1;
        else                                        // (t, j) becomes a mixed pair: add its term
            ppd += w[0] * w[1] - kw[0] * kw[1] + (w[0] - kw[0]) * g0 + (w[1] - kw[1]) * g1;
    }
    const double An0 = Ahat - a[0] * (Aall - a[0]) + exA0, An1 = Ahat - a[1] * (Aall - a[1]) + exA1;
    const double Bn0 = Bhat - b[0] * (Ball - b[0]) + exB0, Bn1 = Bhat - b[1] * (Ball - b[1]) + exB1;
    const double Anp0 = (Aall - a[0]) * (Ball - b[0]) - (Sab - a[0] * b[0]);
    const double Anp1 = (Aall - a[1]) * (Ball - b[1]) - (Sab - a[1] * b[1]);
    const double cs = PA[0] * PB[1] + PA[1] * PB[0] - (KA[0] * KB[1] + KA[1] * KB[0])
                    + (PB[0] - KB[0]) * An0 + (PA[0] - KA[0]) * Bn0 + (PB[1] - KB[1]) * An1 + (PA[1] - KA[1]) * Bn1;
    const double viol = M[0] * M[1] - KM[0] * KM[1] + (M[0] - KM[0]) * Anp0 + (M[1] - KM[1]) * Anp1 - mkboth - eterm;
    const double cnp = cs - viol;
    return (cnp - (pp_[p] + ppd)) - c_[p];
}

// dummy_change at a node p of any degree on j's ancestor path: as the generic arithmetic (branch values
// adjusted for the switch, node_contrib without the dummy-pair terms, j's dummy-pair terms replaced), but
// in one pass over the children and their dummy lists. node_contrib's sums over branches i of products
// with per-branch terms expand into sums over the branches (per dummy t: of x_ti, u_i x_ti, ...) that this
// pass accumulates, so the cost is O(D + sum of the children's list lengths + L) instead of O(D L), and the
// pairs of j with each other dummy only need the branches on j's path.
double WScorer::dummy_change_poly(int p, int j, int tk, int from, double xt, Scratch& S, Scratch::PathVals& me) {
    const SubForest& C = *C_;
    const int lb = dl_off_[p], L = dl_off_[p + 1] - lb, cb = C.child_off[p], D = C.child_off[p + 1] - cb;
    const double sgn = from == 0 ? 1.0 : -1.0;
    // per dummy u of p's list, after the switch: sums over the branches of x, x^2, x b, a x (all branches),
    // x over the child branches; with u_i = PB_i - KB_i, v_i = PA_i - KA_i: of u x, u x^2, v x, v x^2; with
    // wp = Mk - KMk: of wp, wp b, wp a, wp x, wp x b, wp x a; of Mk, Mk^2, KMk, KMk^2 (child branches). For the
    // pairs of j with u, over the child branches holding j: x_u x_j, j's pair sums with u for the parent (W,
    // KW), and the per-pair sums of the pair term (m, m^2, k, k^2, w, w x_j, w x_u, w x_j x_u; m, k: the old
    // mixed-pair sums if (j, u) was mixed, else the new ones)
    enum { gX, gX2, gXb, gXa, gXC, gUx, gUx2, gVx, gVx2, gWP, gWPb, gWPa, gWPx, gWPxb, gWPxa, gM, gM2, gK, gK2,
           gSxj, gWs, gKWs, gPm, gPm2, gPk, gPk2, gPw, gPwj, gPwu, gPwju, gN };
    ensure(S.bx, static_cast<size_t>(gN) * L);
    double* g = S.bx.data();
    std::fill_n(g, static_cast<size_t>(gN) * L, 0.0);
    double* G[gN];
    for (int k = 0; k < gN; ++k) G[k] = g + static_cast<size_t>(k) * L;
    if (S.sd.size() < static_cast<size_t>(L)) S.sd.resize(L);
    int8_t* sd = S.sd.data();                       // sides before the switch
    for (int t = 0; t < L; ++t) sd[t] = part_d_[dl_k_[lb + t]];

    double Aall = 0, Ball = 0, Saa = 0, Sbb = 0, Sab = 0, SPA = 0, SKA = 0, SPB = 0, SKB = 0, SPAPB = 0, SKAKB = 0;
    double SM = 0, SM2 = 0, SKM = 0, SKM2 = 0, U = 0, Ua = 0, Ua2 = 0, V = 0, Vb = 0, Vb2 = 0, Wt = 0, Wa = 0, Wb = 0, Wab = 0;
    double Ssum = 0, Sdot = 0, XjC = 0, Sjj = 0, SSsum = 0, KSSsum = 0;   // j's pairs with its own side (before)
    for (int i = 0; i < D; ++i) {
        const int c = C.children[cb + i];
        const NodeVals& v = nv_[c];
        const double ec = ef_[c], dc = df_[c];
        double a = ec * v.a, b = ec * v.b, PA = v.PA, KA = dc * v.KA, PB = v.PB, KB = dc * v.KB, M = v.M, KM = dc * v.KM;
        const double s_before = from == 0 ? a : b;
        Ssum += s_before;
        const int st = S.ovr[c];                    // -1: j not below c; -3: a cluster of j; >= 0: on j's path
        const double* mkkc = mkk_.data() + pk_off_[c];   // the child's stored pair block
        const double* kmkkc = kmkk_.data() + pk_off_[c];
        const int ob = dl_off_[c], Lc = dl_off_[c + 1] - ob;
        double xj = 0, SSi = 0, KSSi = 0, OPi = 0, KOPi = 0;
        int tjc = -1;
        const Scratch::PathVals* cv = st >= 0 ? &S.pv[st] : nullptr;
        if (st != -1) {
            tjc = static_cast<int>(std::lower_bound(dl_k_.begin() + ob, dl_k_.begin() + ob + Lc, j) - dl_k_.begin()) - ob;
            xj = ec * dl_x_[ob + tjc];
            XjC += xj; Sjj += xj * xj; Sdot += xj * s_before;
            if (from == 0) { a -= xj; b += xj; } else { a += xj; b -= xj; }
            if (cv) {
                SSi = cv->SS; KSSi = dc * cv->KSS; OPi = mk_[ob + tjc]; KOPi = dc * kmk_[ob + tjc];
                SSsum += SSi; KSSsum += KSSi;
                if (from == 0) { PA -= SSi; KA -= KSSi; PB += OPi; KB += KOPi; }
                else { PB -= SSi; KB -= KSSi; PA += OPi; KA += KOPi; }
                M += SSi - OPi;
                KM += KSSi - KOPi;
            }
        }
        const double u = PB - KB, vv = PA - KA, w = M - KM;
        Aall += a; Ball += b; Saa += a * a; Sbb += b * b; Sab += a * b;
        SPA += PA; SKA += KA; SPB += PB; SKB += KB; SPAPB += PA * PB; SKAKB += KA * KB;
        SM += M; SM2 += M * M; SKM += KM; SKM2 += KM * KM;
        U += u; Ua += u * a; Ua2 += u * a * a; V += vv; Vb += vv * b; Vb2 += vv * b * b;
        Wt += w; Wa += w * a; Wb += w * b; Wab += w * a * b;
        int rj = 0, nbc = 0;                          // on j's path: j's rank on its side, the child's B dummies
        if (cv) {
            int na = 0, nb = 0;
            for (int t = 0; t < Lc; ++t) {
                const bool inA = part_d_[dl_k_[ob + t]] == 0;
                if (t == tjc) rj = inA ? na : nb;
                if (inA) ++na; else ++nb;
            }
            nbc = nb;
        }
        for (int t = 0, q = 0, ra = 0, rb = 0; t < Lc; ++t) {
            while (dl_k_[lb + q] != dl_k_[ob + t]) ++q;
            const int r = sd[q] == 0 ? ra++ : rb++;   // its rank among the child's dummies of its side
            const double x = ec * dl_x_[ob + t];
            double mk = mk_[ob + t], kmk = dc * kmk_[ob + t];
            if (cv) {
                if (q == tk) { mk = SSi; kmk = KSSi; }
                else {
                    const double wv = S.W[cv->woff + t], kwv = dc * S.KW[cv->woff + t];
                    mk += (sd[q] == from ? 1.0 : -1.0) * wv;    // same side as j before: (q, j) becomes mixed
                    kmk += (sd[q] == from ? 1.0 : -1.0) * kwv;
                    G[gWs][q] += wv; G[gKWs][q] += kwv;
                    double pm, pk;                              // the pair (j, q): old mixed sums, or new ones
                    if (sd[q] != from) {
                        const size_t sx = from == 0 ? static_cast<size_t>(rj) * nbc + r : static_cast<size_t>(r) * nbc + rj;
                        pm = mkkc[sx]; pk = dc * kmkkc[sx];
                    } else {
                        pm = wv; pk = kwv;
                    }
                    const double pw = pm - pk;
                    G[gPm][q] += pm; G[gPm2][q] += pm * pm; G[gPk][q] += pk; G[gPk2][q] += pk * pk;
                    G[gPw][q] += pw; G[gPwj][q] += pw * xj; G[gPwu][q] += pw * x; G[gPwju][q] += pw * xj * x;
                }
            }
            if (q != tk) G[gSxj][q] += x * xj;
            const double wp = mk - kmk;
            G[gX][q] += x; G[gX2][q] += x * x; G[gXb][q] += x * b; G[gXa][q] += a * x; G[gXC][q] += x;
            G[gUx][q] += u * x; G[gUx2][q] += u * x * x; G[gVx][q] += vv * x; G[gVx2][q] += vv * x * x;
            G[gWP][q] += wp; G[gWPb][q] += wp * b; G[gWPa][q] += wp * a; G[gWPx][q] += wp * x;
            G[gWPxb][q] += wp * x * b; G[gWPxa][q] += wp * x * a;
            G[gM][q] += mk; G[gM2][q] += mk * mk; G[gK][q] += kmk; G[gK2][q] += kmk * kmk;
        }
    }
    // the top branch (no pair sums): j's weight there moves to the other side
    const double aD = aTop_[p] - sgn * xt, bD = bTop_[p] + sgn * xt;
    Aall += aD; Ball += bD; Saa += aD * aD; Sbb += bD * bD; Sab += aD * bD;
    const double* top = &xtop_[xtop_off_[p]];
    for (int t = 0; t < L; ++t) {
        const double x = top[dl_t_[lb + t]];
        G[gX][t] += x; G[gX2][t] += x * x; G[gXb][t] += x * bD; G[gXa][t] += aD * x;
    }

    // j's pair sums at p, for its parent: with each other dummy, and with its own side (before the switch)
    me.woff = static_cast<int>(S.W.size());
    S.W.resize(S.W.size() + L, 0.0);
    S.KW.resize(S.KW.size() + L, 0.0);
    for (int t = 0; t < L; ++t) {
        if (t == tk) continue;
        const double own = G[gXC][t] * XjC - G[gSxj][t];
        S.W[me.woff + t] = own + G[gWs][t];
        S.KW[me.woff + t] = own + G[gKWs][t];
    }
    const double ownS = XjC * Ssum - Sdot - (XjC * XjC - Sjj);
    me.SS = ownS + SSsum;
    me.KSS = ownS + KSSsum;

    // the contribution without the dummy-pair terms, after the switch (node_contrib, expanded)
    double Ahat = (Aall * Aall - Saa) / 2, Bhat = (Ball * Ball - Sbb) / 2;
    double sumA = 0, sumB = 0, dummyM = 0, eterm = 0;
    for (int t = 0; t < L; ++t) {
        const bool inA = (t == tk ? 1 - from : sd[t]) == 0;
        const double X = G[gX][t], own = (X * X - G[gX2][t]) / 2;
        dummyM += (G[gM][t] * G[gM][t] - G[gM2][t]) / 2 - (G[gK][t] * G[gK][t] - G[gK2][t]) / 2;
        if (inA) {
            Ahat -= own;
            sumA += X * G[gUx][t] - G[gUx2][t];
            eterm += (X * Ball - G[gXb][t]) * G[gWP][t] - X * G[gWPb][t] - Ball * G[gWPx][t] + 2 * G[gWPxb][t];
        } else {
            Bhat -= own;
            sumB += X * G[gVx][t] - G[gVx2][t];
            eterm += (Aall * X - G[gXa][t]) * G[gWP][t] - Aall * G[gWPx][t] - X * G[gWPa][t] + 2 * G[gWPxa][t];
        }
    }
    const double cs = SPA * SPB - SPAPB - (SKA * SKB - SKAKB)
                    + Ahat * U - Aall * Ua + Ua2 + sumA + Bhat * V - Ball * Vb + Vb2 + sumB;
    const double viol = (SM * SM - SM2) / 2 - (SKM * SKM - SKM2) / 2 - dummyM
                      + Wt * (Aall * Ball - Sab) - Ball * Wa - Aall * Wb + 2 * Wab - eterm;
    const double cnp = cs - viol;

    // dummy-pair terms: the stored ones, minus j's old pairs, plus its new ones (pair_term from the sums)
    double pp = pp_[p];
    const double Xj = G[gX][tk];
    for (int t = 0; t < L; ++t) {
        if (t == tk || (G[gPw][t] == 0 && G[gPm][t] == 0 && G[gPk][t] == 0)) continue;
        const double Xu = G[gX][t], sxy = G[gSxj][t] + top[dl_t_[lb + t]] * xt;
        const double term = (G[gPm][t] * G[gPm][t] - G[gPm2][t]) / 2 - (G[gPk][t] * G[gPk][t] - G[gPk2][t]) / 2
                          + G[gPw][t] * (Xj * Xu - sxy) - Xj * G[gPwu][t] - Xu * G[gPwj][t] + 2 * G[gPwju][t];
        pp += sd[t] != from ? -term : term;
    }
    return (cnp - pp) - c_[p];
}

void WScorer::set_partition(const std::vector<int8_t>& part_r, const std::vector<int8_t>& part_d) {
    part_r_ = part_r;
    part_d_ = part_d;
}

double WScorer::score() {
    run_chunks([&](int c, int slot) {
        double s = 0;
        for (int g = chunk_[c]; g < chunk_[c + 1]; ++g) s += tree_pass(g, part_r_.data(), part_d_.data(), false, scratch_[slot]);
        chunk_score_[c] = s;
    });
    double s = 0;
    for (double x : chunk_score_) s += x;
    return s;
}

void WScorer::scores(const ScoreBatch& batch, std::vector<double>& out) {
    const int Q = batch.size(), Cn = static_cast<int>(chunk_.size()) - 1, W = R_ + d_;
    out.assign(Q, 0.0);
    if (Q == 0) return;
    batch_score_.resize(static_cast<size_t>(Q) * Cn);
    run_chunks([&](int c, int slot) {
        for (int q = 0; q < Q; ++q) {
            const int8_t* pr = batch.sides.data() + static_cast<size_t>(q) * W;
            double s = 0;
            for (int g = chunk_[c]; g < chunk_[c + 1]; ++g) s += tree_pass(g, pr, pr + R_, false, scratch_[slot]);
            batch_score_[static_cast<size_t>(q) * Cn + c] = s;
        }
    });
    const int8_t* last = batch.sides.data() + static_cast<size_t>(Q - 1) * W;
    part_r_.assign(last, last + R_);
    part_d_.assign(last + R_, last + W);
    for (int q = 0; q < Q; ++q) {
        double s = 0;
        for (int c = 0; c < Cn; ++c) s += batch_score_[static_cast<size_t>(q) * Cn + c];
        out[q] = s;
    }
}

double WScorer::gains(std::vector<double>& real_gain_out, std::vector<double>& dummy_gain) {
    const int W = R_ + d_;
    run_chunks([&](int c, int slot) {
        Scratch& S = scratch_[slot];
        double* rg = &chunk_gain_[static_cast<size_t>(c) * W];
        std::fill(rg, rg + W, 0.0);
        double s = 0;
        for (int g = chunk_[c]; g < chunk_[c + 1]; ++g) {
            const double sg = tree_pass(g, part_r_.data(), part_d_.data(), true, S);
            s += sg;
            gain_pass(g, part_r_.data(), part_d_.data(), rg, S);
            if (pres_off_[g + 1] > pres_off_[g]) {
                S.pd = part_d_;
                for (int q = pres_off_[g]; q < pres_off_[g + 1]; ++q)
                    rg[R_ + pres_[q]] += dummy_change(g, q - pres_off_[g], pres_[q], S);
            }
        }
        chunk_score_[c] = s;
    });
    real_gain_out.assign(R_, 0.0);
    dummy_gain.assign(d_, 0.0);
    double s = 0;
    for (size_t c = 0; c < chunk_score_.size(); ++c) {
        s += chunk_score_[c];
        const double* rg = &chunk_gain_[c * W];
        for (int i = 0; i < R_; ++i) real_gain_out[i] += rg[i];
        for (int j = 0; j < d_; ++j) dummy_gain[j] += rg[R_ + j];
    }
    return s;
}

}  // namespace wqfm
