// Brute-force validation of the Scorer.
//
// Score(A,B,G) = sum_g sum_q w(q) * (+1 satisfied, -1 violated, 0 unresolved), over quartets with two
// taxa on each side and no two taxa under the same dummy (paper Section 2.5). Random gene trees
// (polytomies, missing taxa, rooted and unrooted) and random subproblems with nested dummy taxa;
// checks the score, every real and dummy gain, and the state after incremental swaps.
//
// usage: brute_check [seed] [trials]     exits non-zero on any mismatch

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "forest.hpp"
#include "newick.hpp"
#include "scorer.hpp"
#include "units.hpp"

using namespace wqfm;

namespace {

std::mt19937_64 rng;

int rnd(int n) { return static_cast<int>(rng() % static_cast<uint64_t>(n)); }
double unif() { return std::uniform_real_distribution<double>(0, 1)(rng); }

std::string random_tree(std::vector<std::string> parts, double poly, bool rooted) {
    std::shuffle(parts.begin(), parts.end(), rng);
    const size_t stop = rooted ? 2 : 3;
    while (parts.size() > stop) {
        int k = 2;
        if (unif() < poly && parts.size() > stop + 2) k = 3 + rnd(std::min<int>(3, static_cast<int>(parts.size() - stop - 1)));
        std::string node = "(";
        for (int i = 0; i < k; ++i) {
            const int at = rnd(static_cast<int>(parts.size()));
            node += (i ? "," : "") + parts[at];
            parts.erase(parts.begin() + at);
        }
        parts.push_back(node + ")");
    }
    std::string s = "(";
    for (size_t i = 0; i < parts.size(); ++i) s += (i ? "," : "") + parts[i];
    return s + ");";
}

struct Brute {
    const Forest& F;
    std::vector<std::vector<int>> dist;   // per tree: leaf-to-leaf distances, indexed [taxon*n + taxon]

    explicit Brute(const Forest& f) : F(f) {
        for (int g = 0; g < F.k; ++g) {
            const int off = F.tree_off[g], m = F.tree_off[g + 1] - off;
            std::vector<std::vector<int>> adj(m);
            for (int v = 0; v < m; ++v)
                if (F.parent[off + v] >= 0) { adj[v].push_back(F.parent[off + v] - off); adj[F.parent[off + v] - off].push_back(v); }
            std::vector<int> d(static_cast<size_t>(F.n) * F.n, -1);
            for (int v = 0; v < m; ++v) {
                if (F.taxon[off + v] < 0) continue;
                std::vector<int> dv(m, -1);
                std::queue<int> q;
                dv[v] = 0; q.push(v);
                while (!q.empty()) {
                    const int u = q.front(); q.pop();
                    for (int w : adj[u]) if (dv[w] < 0) { dv[w] = dv[u] + 1; q.push(w); }
                }
                for (int u = 0; u < m; ++u)
                    if (F.taxon[off + u] >= 0) d[static_cast<size_t>(F.taxon[off + v]) * F.n + F.taxon[off + u]] = dv[u];
            }
            dist.push_back(std::move(d));
        }
    }

    // side[t], weight[t], dummy[t] (-1 for real taxa)
    double score(const std::vector<int>& side, const std::vector<double>& weight, const std::vector<int>& dummy) const {
        double total = 0;
        auto same = [&](int x, int y) { return dummy[x] >= 0 && dummy[x] == dummy[y]; };
        for (int g = 0; g < F.k; ++g) {
            const std::vector<int>& d = dist[g];
            auto D = [&](int x, int y) { return d[static_cast<size_t>(x) * F.n + y]; };
            std::vector<int> A, B;
            for (int t = 0; t < F.n; ++t) if (F.leaf(g, t) >= 0) (side[t] == 0 ? A : B).push_back(t);
            for (size_t i = 0; i < A.size(); ++i)
                for (size_t j = i + 1; j < A.size(); ++j) {
                    const int a = A[i], b = A[j];
                    if (same(a, b)) continue;
                    for (size_t k = 0; k < B.size(); ++k)
                        for (size_t l = k + 1; l < B.size(); ++l) {
                            const int c = B[k], e = B[l];
                            if (same(c, e) || same(a, c) || same(a, e) || same(b, c) || same(b, e)) continue;
                            const int s1 = D(a, b) + D(c, e), s2 = D(a, c) + D(b, e), s3 = D(a, e) + D(b, c);
                            const int topo = (s1 < s2 && s1 < s3) ? 1 : (s1 == s2 && s2 == s3) ? 0 : -1;
                            total += weight[a] * weight[b] * weight[c] * weight[e] * topo;
                        }
                }
        }
        return total;
    }
};

double rel(double x, double y) { return std::fabs(x - y) / std::max(1.0, std::fabs(y)); }

}  // namespace

int main(int argc, char** argv) {
    const uint64_t seed = argc > 1 ? std::stoull(argv[1]) : 1;
    const int trials = argc > 2 ? std::stoi(argv[2]) : 300;
    rng.seed(seed);
    double worst_score = 0, worst_gain = 0;
    int checks = 0;

    for (int trial = 0; trial < trials; ++trial) {
        const int n = 6 + rnd(10), k = 1 + rnd(6);
        const double poly = unif() < 0.5 ? 0.0 : 0.4;
        Taxa taxa;
        std::vector<std::string> labels;
        for (int i = 0; i < n; ++i) { labels.push_back("t" + std::to_string(i)); taxa.add(labels.back()); }
        std::vector<ParsedTree> trees;
        for (int g = 0; g < k; ++g) {
            std::vector<std::string> present = labels;
            if (g > 0 && unif() < 0.4) {
                const int drop = rnd(std::max(1, n - 5));
                for (int i = 0; i < drop; ++i) present.erase(present.begin() + rnd(static_cast<int>(present.size())));
            }
            trees.push_back(parse_newick(random_tree(present, poly, unif() < 0.3)));
        }
        const Forest F = make_forest(trees, taxa);
        const Brute brute(F);

        // random subproblem with nested dummies
        DummyPool pool;
        std::vector<int> order(n);
        for (int i = 0; i < n; ++i) order[i] = i;
        std::shuffle(order.begin(), order.end(), rng);
        std::vector<int> dummies;
        const int nd = rnd(std::min(4, n / 3) + 1);
        const bool all_dummies = rnd(5) == 0;       // deep subproblems can have no real taxa at all
        for (int j = 0; j < nd; ++j) {
            std::vector<int> rs;
            const int take = 1 + rnd(3);
            for (int i = 0; i < take && order.size() > 4; ++i) { rs.push_back(order.back()); order.pop_back(); }
            if (rs.empty()) break;
            std::vector<int> subs;
            if (!dummies.empty() && rnd(2)) { const int at = rnd(static_cast<int>(dummies.size())); subs.push_back(dummies[at]); dummies.erase(dummies.begin() + at); }
            dummies.push_back(pool.add(rs, subs));
        }
        if (all_dummies) {                          // remaining taxa become dummies of 1-3 taxa each
            while (!order.empty()) {
                std::vector<int> rs;
                const int take = 1 + rnd(3);
                for (int i = 0; i < take && !order.empty(); ++i) { rs.push_back(order.back()); order.pop_back(); }
                dummies.push_back(pool.add(rs, {}));
            }
        }
        Subproblem sp;
        sp.reals = order;
        sp.dummies = dummies;
        if (sp.size() < 2) continue;
        for (size_t i = 0; i < sp.reals.size(); ++i) sp.part_r.push_back(static_cast<int8_t>(rnd(2)));
        for (size_t j = 0; j < sp.dummies.size(); ++j) sp.part_d.push_back(static_cast<int8_t>(rnd(2)));

        // independent weights: product of child counts on the path (Eq. 1)
        std::vector<double> weight(n, 1.0);
        std::vector<int> dummy_of(n, -1);
        auto assign = [&](auto&& self, int id, double denom, int top) -> void {
            const Dummy& d = pool[id];
            const double c = denom * static_cast<double>(d.reals.size() + d.subs.size());
            for (int r : d.reals) { weight[r] = 1.0 / c; dummy_of[r] = top; }
            for (int s : d.subs) self(self, s, c, top);
        };
        for (size_t j = 0; j < sp.dummies.size(); ++j) assign(assign, sp.dummies[j], 1.0, static_cast<int>(j));
        auto sides = [&](const std::vector<int8_t>& pr, const std::vector<int8_t>& pd) {
            std::vector<int> side(n);
            for (size_t i = 0; i < sp.reals.size(); ++i) side[sp.reals[i]] = pr[i];
            for (int t = 0; t < n; ++t) if (dummy_of[t] >= 0) side[t] = pd[dummy_of[t]];
            return side;
        };

        Scorer sc(F);
        sc.setup(sp, pool);
        sc.set_partition(sp.part_r, sp.part_d);
        std::vector<double> gr, gd;
        const double fast = sc.gains(gr, gd);
        const double base = brute.score(sides(sp.part_r, sp.part_d), weight, dummy_of);
        worst_score = std::max(worst_score, rel(fast, base));
        worst_score = std::max(worst_score, rel(sc.score(), base));
        for (size_t i = 0; i < sp.reals.size(); ++i) {
            auto pr = sp.part_r; pr[i] ^= 1;
            const double expect = brute.score(sides(pr, sp.part_d), weight, dummy_of) - base;
            const double e = rel(gr[i], expect);
            worst_gain = std::max(worst_gain, e);
            if (e > 1e-9) std::printf("REAL GAIN trial %d unit %zu: %.10g vs brute %.10g\n", trial, i, gr[i], expect);
        }
        for (size_t j = 0; j < sp.dummies.size(); ++j) {
            auto pd = sp.part_d; pd[j] ^= 1;
            const double expect = brute.score(sides(sp.part_r, pd), weight, dummy_of) - base;
            const double e = rel(gd[j], expect);
            worst_gain = std::max(worst_gain, e);
            if (e > 1e-9) std::printf("DUMMY GAIN trial %d dummy %zu: %.10g vs brute %.10g\n", trial, j, gd[j], expect);
        }
        // incremental swaps, then score and gains again
        auto pr = sp.part_r, pd = sp.part_d;
        if (!pr.empty()) { const int i = rnd(static_cast<int>(pr.size())); sc.swap_real(i); pr[i] ^= 1; }
        if (!pd.empty()) { const int j = rnd(static_cast<int>(pd.size())); sc.swap_dummy(j); pd[j] ^= 1; }
        const double after = brute.score(sides(pr, pd), weight, dummy_of);
        const double fast2 = sc.gains(gr, gd);
        worst_score = std::max(worst_score, rel(fast2, after));
        for (size_t i = 0; i < pr.size(); ++i) {
            auto p2 = pr; p2[i] ^= 1;
            worst_gain = std::max(worst_gain, rel(gr[i], brute.score(sides(p2, pd), weight, dummy_of) - after));
        }
        for (size_t j = 0; j < pd.size(); ++j) {
            auto p2 = pd; p2[j] ^= 1;
            worst_gain = std::max(worst_gain, rel(gd[j], brute.score(sides(pr, p2), weight, dummy_of) - after));
        }
        ++checks;
    }
    std::printf("brute_check seed=%llu subproblems=%d  max rel err: score %.2e, gains %.2e\n",
                static_cast<unsigned long long>(seed), checks, worst_score, worst_gain);
    return (worst_score > 1e-9 || worst_gain > 1e-9) ? 1 : 0;
}
