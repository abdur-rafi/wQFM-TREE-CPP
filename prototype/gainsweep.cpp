// Prototype of the wQFM-TREE top-level gain sweep (no dummy taxa), for benchmarking against the Java code.
//
// Per node u with branches (children + top), a_i / b_i = number of partition-0 / partition-1 taxa in branch i.
// Node contribution f_u = w(S(g,u)) + w(U(g,u))/2 (the Java "sat" per node), score = sum_g (2 sum_u f_u - w(PA(g)) w(PB(g))).
// Binary nodes use a closed form; polytomies use power sums so all branch gains cost O(D).
// Gains for every taxon come from one bottom-up pass (per-node deltas) and one top-down pass (path sums).
// At the top level every quantity is an integer < 2^53, so doubles are exact and results match Java bit for bit.

#include <algorithm>
#include <atomic>
#include <barrier>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using i32 = int32_t;

struct Forest {
    int n = 0, k = 0;
    std::vector<i32> tree_off;   // k+1
    std::vector<i32> parent;     // -1 for a tree root
    std::vector<i32> child_off;  // N+1, CSR into children
    std::vector<i32> children;
    std::vector<i32> taxon;      // -1 for internal nodes
    std::vector<i32> leaf_of;    // k*n, -1 if taxon absent from tree
};

static Forest read_forest(const std::string& path, const std::unordered_map<std::string, int>& ids) {
    Forest F;
    F.n = (int)ids.size();
    std::ifstream in(path);
    std::string line;
    F.tree_off.push_back(0);
    F.child_off.push_back(0);
    while (std::getline(in, line)) {
        if (line.find('(') == std::string::npos) continue;
        std::vector<std::vector<i32>> stack;
        size_t i = 0, L = line.size();
        auto new_node = [&](i32 tax, const std::vector<i32>* kids) {
            i32 v = (i32)F.parent.size();
            F.parent.push_back(-1);
            F.taxon.push_back(tax);
            if (kids) for (i32 c : *kids) { F.children.push_back(c); F.parent[c] = v; }
            F.child_off.push_back((i32)F.children.size());
            return v;
        };
        auto skip_suffix = [&]() {  // internal label / support, then optional :length
            while (i < L && line[i] != ',' && line[i] != ')' && line[i] != ';' && line[i] != ':' && line[i] != '(') ++i;
            if (i < L && line[i] == ':') { ++i; while (i < L && line[i] != ',' && line[i] != ')' && line[i] != ';') ++i; }
        };
        i32 root = -1;
        while (i < L) {
            char c = line[i];
            if (c == '(') { stack.emplace_back(); ++i; }
            else if (c == ')') {
                ++i;
                std::vector<i32> kids = std::move(stack.back());
                stack.pop_back();
                i32 v = new_node(-1, &kids);
                skip_suffix();
                if (stack.empty()) root = v; else stack.back().push_back(v);
            } else if (c == ',' || c == ';' || isspace((unsigned char)c)) { ++i; }
            else {
                size_t j = i;
                while (j < L && line[j] != ',' && line[j] != ')' && line[j] != ':' && line[j] != ';') ++j;
                std::string lab = line.substr(i, j - i);
                auto it = ids.find(lab);
                if (it == ids.end()) { std::fprintf(stderr, "unknown taxon %s\n", lab.c_str()); std::exit(1); }
                i = j;
                i32 v = new_node(it->second, nullptr);
                skip_suffix();
                stack.back().push_back(v);
            }
        }
        if (root != (i32)F.parent.size() - 1) { std::fprintf(stderr, "parse error\n"); std::exit(1); }
        F.tree_off.push_back((i32)F.parent.size());
        F.k++;
    }
    F.leaf_of.assign((size_t)F.k * F.n, -1);
    for (int g = 0; g < F.k; ++g)
        for (i32 v = F.tree_off[g]; v < F.tree_off[g + 1]; ++v)
            if (F.taxon[v] >= 0) F.leaf_of[(size_t)g * F.n + F.taxon[v]] = v;
    return F;
}

// ---------- node formulas ----------
struct Agg { double A, B, Saa, Sbb, Sp, Sap, Saap, Sab, Sabb, Saab, Saabb; };

static inline void agg_add(Agg& g, double a, double b, double s) {
    double p = b * (b - 1) / 2;
    g.A += s * a; g.B += s * b; g.Saa += s * a * a; g.Sbb += s * b * b; g.Sp += s * p;
    g.Sap += s * a * p; g.Saap += s * a * a * p; g.Sab += s * a * b; g.Sabb += s * a * b * b;
    g.Saab += s * a * a * b; g.Saabb += s * a * a * b * b;
}

// f = w(S) + w(U)/2 from power sums (top level: PA_ij = a_i a_j, PB_k = C(b_k,2), PB_kl = b_k b_l)
static inline double agg_f(const Agg& g) {
    double E2a = (g.A * g.A - g.Saa) / 2, E2b = (g.B * g.B - g.Sbb) / 2;
    double S = g.Sp * E2a - g.A * g.Sap + g.Saap;
    double U = E2a * E2b - (g.A * (g.B * g.Sab - g.Sabb) - (g.B * g.Saab - g.Saabb)) + (g.Sab * g.Sab - g.Saabb) / 2;
    return S + U / 2;
}

struct State {
    const Forest* F;
    std::vector<i32> part;              // taxon -> 0/1
    std::vector<i32> below[2];          // per node
    std::vector<i32> tot[2];            // per tree
    std::vector<double> acc;            // per node, 2 entries: push values, then path sums
    int cnt[2] = {0, 0};

    void init(const Forest& f, const std::vector<i32>& p) {
        F = &f; part = p;
        size_t N = f.parent.size();
        below[0].assign(N, 0); below[1].assign(N, 0);
        tot[0].assign(f.k, 0); tot[1].assign(f.k, 0);
        acc.assign(2 * N, 0.0);
        cnt[0] = cnt[1] = 0;
        for (int t = 0; t < f.n; ++t) cnt[part[t]]++;
        for (int g = 0; g < f.k; ++g) {
            for (i32 v = f.tree_off[g]; v < f.tree_off[g + 1]; ++v) {
                if (f.taxon[v] >= 0) below[part[f.taxon[v]]][v] = 1;
                else for (i32 e = f.child_off[v]; e < f.child_off[v + 1]; ++e) {
                    i32 c = f.children[e];
                    below[0][v] += below[0][c]; below[1][v] += below[1][c];
                }
            }
            i32 root = f.tree_off[g + 1] - 1;
            tot[0][g] = below[0][root]; tot[1][g] = below[1][root];
        }
    }

    // processes trees [g0,g1); adds gains (own-partition direction) into gain[taxon]; returns score part
    double sweep(int g0, int g1, double* gain) {
        const Forest& f = *F;
        const i32* b0 = below[0].data(); const i32* b1 = below[1].data();
        std::vector<double> da, db;
        double score = 0;
        for (int g = g0; g < g1; ++g) {
            const double T0 = tot[0][g], T1 = tot[1][g];
            double G0 = 0, G1 = 0, fsum = 0;
            const i32 off = f.tree_off[g], end = f.tree_off[g + 1];
            for (i32 v = off; v < end; ++v) {
                if (f.taxon[v] >= 0) continue;
                const i32 cb = f.child_off[v], ce = f.child_off[v + 1], nc = ce - cb;
                const double at = T0 - b0[v], bt = T1 - b1[v];
                if (nc == 2) {
                    const i32 c0 = f.children[cb], c1 = f.children[cb + 1];
                    const double a[3] = {(double)b0[c0], (double)b0[c1], at};
                    const double b[3] = {(double)b1[c0], (double)b1[c1], bt};
                    const double p[3] = {b[0] * (b[0] - 1) / 2, b[1] * (b[1] - 1) / 2, b[2] * (b[2] - 1) / 2};
                    fsum += a[0] * a[1] * p[2] + a[1] * a[2] * p[0] + a[0] * a[2] * p[1];
                    double dAB[3], dBA[3];
                    for (int i = 0; i < 3; ++i) {
                        const int j = (i + 1) % 3, l = (i + 2) % 3;
                        const double cross = a[j] * p[l] + a[l] * p[j];
                        dAB[i] = a[i] > 0 ? a[j] * a[l] * b[i] - cross : 0.0;
                        dBA[i] = b[i] > 0 ? cross - a[j] * a[l] * (b[i] - 1) : 0.0;
                    }
                    acc[2 * c0] = dAB[0] - dAB[2]; acc[2 * c0 + 1] = dBA[0] - dBA[2];
                    acc[2 * c1] = dAB[1] - dAB[2]; acc[2 * c1 + 1] = dBA[1] - dBA[2];
                    G0 += dAB[2]; G1 += dBA[2];
                } else {
                    Agg A{};
                    for (i32 e = cb; e < ce; ++e) { i32 c = f.children[e]; agg_add(A, b0[c], b1[c], 1); }
                    agg_add(A, at, bt, 1);
                    const double f0 = agg_f(A);
                    fsum += f0;
                    auto delta = [&](double a, double b, int dir) {
                        Agg h = A;
                        agg_add(h, a, b, -1);
                        if (dir == 0) agg_add(h, a - 1, b + 1, 1); else agg_add(h, a + 1, b - 1, 1);
                        return agg_f(h) - f0;
                    };
                    const double tAB = at > 0 ? delta(at, bt, 0) : 0.0, tBA = bt > 0 ? delta(at, bt, 1) : 0.0;
                    for (i32 e = cb; e < ce; ++e) {
                        i32 c = f.children[e];
                        double x0 = b0[c], x1 = b1[c];
                        acc[2 * c] = (x0 > 0 ? delta(x0, x1, 0) : 0.0) - tAB;
                        acc[2 * c + 1] = (x1 > 0 ? delta(x0, x1, 1) : 0.0) - tBA;
                    }
                    G0 += tAB; G1 += tBA;
                }
            }
            // top-down path sums (parents have larger indices)
            const i32 root = end - 1;
            acc[2 * root] = 0; acc[2 * root + 1] = 0;
            const double PA = T0 * (T0 - 1) / 2, PB = T1 * (T1 - 1) / 2, tq = PA * PB;
            const double dTot0 = ((T0 - 1) * (T0 - 2) / 2) * ((T1 + 1) * T1 / 2) - tq;
            const double dTot1 = ((T0 + 1) * T0 / 2) * ((T1 - 1) * (T1 - 2) / 2) - tq;
            for (i32 v = root - 1; v >= off; --v) {
                const i32 pv = f.parent[v];
                acc[2 * v] += acc[2 * pv]; acc[2 * v + 1] += acc[2 * pv + 1];
                const i32 t = f.taxon[v];
                if (t >= 0) {
                    if (part[t] == 0) gain[t] += 2 * (acc[2 * v] + G0) - dTot0;
                    else gain[t] += 2 * (acc[2 * v + 1] + G1) - dTot1;
                }
            }
            score += 2 * fsum - tq;
        }
        return score;
    }

    void swap_taxon_trees(int y, int g0, int g1) {
        const Forest& f = *F;
        const int p = part[y], q = 1 - p;
        for (int g = g0; g < g1; ++g) {
            i32 v = f.leaf_of[(size_t)g * f.n + y];
            if (v < 0) continue;
            for (; v >= 0; v = f.parent[v]) { below[p][v]--; below[q][v]++; }
            tot[p][g]--; tot[q][g]++;
        }
    }
};

// ---------- tiny persistent fork-join pool ----------
struct Pool {
    int T;
    std::vector<std::thread> th;
    std::barrier<> start_b, done_b;
    std::function<void(int)> job;
    bool stop = false;
    explicit Pool(int t) : T(t), start_b(t), done_b(t) {
        for (int i = 1; i < T; ++i)
            th.emplace_back([this, i] {
                for (;;) { start_b.arrive_and_wait(); if (stop) return; job(i); done_b.arrive_and_wait(); }
            });
    }
    void run(std::function<void(int)> f) {
        job = std::move(f);
        if (T > 1) start_b.arrive_and_wait();
        job(0);
        if (T > 1) done_b.arrive_and_wait();
    }
    ~Pool() { stop = true; if (T > 1) start_b.arrive_and_wait(); for (auto& x : th) x.join(); }
};

struct Engine {
    const Forest& F; State S; Pool pool; int T;
    // fixed chunking (independent of T) + per-chunk buffers summed in chunk order -> identical results for any T;
    // chunks are claimed dynamically so slow (efficiency) cores and polytomy-heavy trees do not stall the sweep
    std::vector<int> chunk;                 // C+1 tree boundaries, balanced by node count
    std::vector<std::vector<double>> cgain;
    std::vector<double> cscore;
    std::atomic<int> next{0};
    Engine(const Forest& f, const std::vector<i32>& part, int t, int C = 64) : F(f), pool(t), T(t) {
        S.init(f, part);
        C = std::max(1, std::min(C, f.k));
        chunk.assign(C + 1, f.k);
        chunk[0] = 0;
        size_t N = f.parent.size();
        for (int b = 1; b < C; ++b) {
            size_t target = N * b / C;
            chunk[b] = (int)(std::lower_bound(f.tree_off.begin(), f.tree_off.end(), (i32)target) - f.tree_off.begin());
            chunk[b] = std::max(chunk[b], chunk[b - 1]);
        }
        cgain.assign(C, std::vector<double>(f.n));
        cscore.assign(C, 0);
    }
    double gains(std::vector<double>& out) {
        const int C = (int)cgain.size();
        next.store(0);
        pool.run([&](int) {
            for (int c; (c = next.fetch_add(1)) < C;) {
                std::fill(cgain[c].begin(), cgain[c].end(), 0.0);
                cscore[c] = S.sweep(chunk[c], chunk[c + 1], cgain[c].data());
            }
        });
        out.assign(F.n, 0.0);
        double sc = 0;
        for (int c = 0; c < C; ++c) {
            sc += cscore[c];
            const double* gp = cgain[c].data();
            for (int t = 0; t < F.n; ++t) out[t] += gp[t];
        }
        return sc;
    }
    void swap(int y) {                      // O(depth) per tree: cheaper serially than a fork-join round trip
        S.swap_taxon_trees(y, 0, F.k);
        S.cnt[S.part[y]]--; S.part[y] = 1 - S.part[y]; S.cnt[S.part[y]]++;
    }
};

int main(int argc, char** argv) {
    if (argc < 8) { std::fprintf(stderr, "usage: gainsweep genetrees taxa_order partition threads reps mode(sweep|fm) outprefix\n"); return 1; }
    std::string gt = argv[1], taxf = argv[2], partf = argv[3], mode = argv[6], outp = argv[7];
    int T = std::atoi(argv[4]), reps = std::atoi(argv[5]);

    std::unordered_map<std::string, int> ids;
    { std::ifstream in(taxf); int id; std::string lab; while (in >> id >> lab) ids[lab] = id; }
    std::vector<i32> part(ids.size());
    { std::ifstream in(partf); int id, p; while (in >> id >> p) part[id] = p; }

    auto t0 = std::chrono::steady_clock::now();
    Forest F = read_forest(gt, ids);
    auto t1 = std::chrono::steady_clock::now();
    size_t N = F.parent.size();
    std::fprintf(stderr, "trees=%d taxa=%d nodes=%zu parse=%.1f ms threads=%d\n", F.k, F.n, N,
                 std::chrono::duration<double, std::milli>(t1 - t0).count(), T);

    Engine E(F, part, T);
    std::vector<double> g;
    if (mode == "sweep") {
        double sc = E.gains(g);            // warm-up
        auto s0 = std::chrono::steady_clock::now();
        for (int r = 0; r < reps; ++r) sc = E.gains(g);
        auto s1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(s1 - s0).count() / reps;
        std::printf("score %.17g\nsweep_ms %.4f\nns_per_node %.3f\n", sc, ms, ms * 1e6 / N);
        FILE* fo = std::fopen((outp + ".cpp.gains").c_str(), "w");
        for (int t = 0; t < F.n; ++t) std::fprintf(fo, "%d %.17g\n", t, g[t]);
        std::fclose(fo);
    } else {
        // one FM pass, replicating QFM2.oneInteration / swapMax at the top level (no dummies, singletons disallowed)
        auto s0 = std::chrono::steady_clock::now();
        const double EPS = 1e-5;
        double cg = 0, maxCg = 0; int maxIdx = -1;
        bool singleton = E.S.cnt[0] == 1 || E.S.cnt[1] == 1;
        std::vector<char> locked(F.n, 0);
        std::vector<int> swaps;
        FILE* fo = std::fopen((outp + ".cpp.fmlog").c_str(), "w");
        for (;;) {
            E.gains(g);
            int best = -1; double bg = 0;
            for (int t = 0; t < F.n; ++t) {
                if (locked[t]) continue;
                if (E.S.cnt[E.S.part[t]] > 2 && (best == -1 || bg < g[t])) { best = t; bg = g[t]; }
            }
            if (best == -1) break;
            E.swap(best); locked[best] = 1; swaps.push_back(best);
            cg += bg;
            std::fprintf(fo, "%d %.17g %.17g\n", best, bg, cg);
            if (singleton && maxIdx == -1) { maxCg = cg; maxIdx = (int)swaps.size() - 1; }
            if (cg > maxCg && std::fabs(maxCg - cg) > EPS) { maxCg = cg; maxIdx = (int)swaps.size() - 1; }
        }
        bool improved = maxIdx != -1;
        int keep = improved ? maxIdx + 1 : ((int)swaps.size() == F.n ? (int)swaps.size() : 0);
        for (int i = (int)swaps.size() - 1; i >= keep; --i) E.swap(swaps[i]);
        auto s1 = std::chrono::steady_clock::now();
        std::fclose(fo);
        double sc = E.gains(g);
        std::printf("fm_moves %zu kept %d improved %d max_cg %.17g final_score %.17g\nfm_pass_s %.3f\n",
                    swaps.size(), keep, (int)improved, maxCg, sc, std::chrono::duration<double>(s1 - s0).count());
        FILE* fp = std::fopen((outp + ".cpp.fmpart").c_str(), "w");
        for (int t = 0; t < F.n; ++t) std::fprintf(fp, "%d %d\n", t, E.S.part[t]);
        std::fclose(fp);
    }
    return 0;
}
