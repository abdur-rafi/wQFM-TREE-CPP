// Brute-force validation of the Scorer.
//
// Score(A,B,G) = sum_g sum_q w(q) * (+1 satisfied, -1 violated, 0 unresolved), over quartets with two
// taxa on each side and no two taxa under the same dummy (paper Section 2.5). Random gene trees
// (polytomies, missing taxa, rooted and unrooted) and random subproblems with nested dummy taxa;
// checks the score, every real and dummy gain, and the state after incremental swaps, for both the
// compressed and the plain subproblem forest. Also checks resolve_polytomies (binary trees on the same
// taxa that keep every resolved quartet) and quartet_support (exact counts, vs enumeration).
//
// usage: brute_check [seed] [trials] [kappa] [gamma] [weighting]     exits non-zero on any mismatch
//   kappa, gamma: the experimental objective variants (Scorer::set_kappa, set_dummy_gamma); a quartet
//   then counts 2 - kappa if satisfied, 1 - kappa if unresolved and -kappa if violated
//   weighting: 0 none, 1 support, 2 length, 3 both. The reference weights a resolved quartet as in
//   weighted ASTRAL: 1 - prod(1 - s_e) over the branches e between its two anchor nodes (a branch
//   without support counts as 1), times exp(-L), L the length of the four paths from the leaves to
//   their anchors (a branch without length counts as 0). In these modes the weighted scorer (WScorer)
//   is checked against the reference, and the reference against itself.

#include <algorithm>
#include <cmath>
#include <memory>
#include <cstdio>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "forest.hpp"
#include "newick.hpp"
#include "qscore.hpp"
#include "scorer.hpp"
#include "wscorer.hpp"
#include "units.hpp"

using namespace wqfm;

namespace {

std::mt19937_64 rng;
double g_kappa = 1, g_gamma = 0;
double dummy_weight(double c) { return g_gamma == 0 ? 1.0 / c : std::pow(c, g_gamma - 1); }

int rnd(int n) { return static_cast<int>(rng() % static_cast<uint64_t>(n)); }
double unif() { return std::uniform_real_distribution<double>(0, 1)(rng); }
int g_weight = 0;
std::string num(double x) { char b[32]; std::snprintf(b, sizeof b, "%.4g", x); return b; }
std::string random_length() { return unif() < 0.1 ? "" : ":" + num(unif() * 1.5); }
std::string random_support() {                 // "" (none), exactly 0 or 1 often, else uniform
    const double u = unif();
    return u < 0.1 ? "" : u < 0.25 ? "1" : u < 0.35 ? "0" : num(unif());
}

// Newick with random supports and branch lengths (some absent)
std::string random_tree(std::vector<std::string> parts, double poly, bool rooted) {
    std::shuffle(parts.begin(), parts.end(), rng);
    for (auto& p : parts) p += random_length();
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
        parts.push_back(node + ")" + random_support() + random_length());
    }
    std::string s = "(";
    for (size_t i = 0; i < parts.size(); ++i) s += (i ? "," : "") + parts[i];
    return s + ");";
}

struct Brute {
    const Forest& F;
    int mode;                             // weighting: 0 none, 1 support, 2 length, 3 both
    std::vector<std::vector<int>> dist;   // per tree: leaf-to-leaf distances, indexed [taxon*n + taxon]
    struct Edge { double s, l; };
    std::vector<std::vector<Edge>> edges;      // per tree: the unrooted tree's edges (a degree-2 root removed)
    std::vector<std::vector<uint64_t>> pmask;  // per tree: edges on the path between two taxa, [taxon*n + taxon]

    explicit Brute(const Forest& f, int weighting = 0) : F(f), mode(weighting) {
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

            // the unrooted tree: every node's branch to its parent is an edge, except that the two
            // branches of a degree-2 root are one edge
            std::vector<std::vector<std::pair<int, int>>> ua(m);   // (neighbour, edge)
            std::vector<Edge> es;
            const int root = m - 1;
            const bool root2 = F.child_off[off + root + 1] - F.child_off[off + root] == 2;
            for (int v = 0; v < root; ++v) {
                const int p = F.parent[off + v] - off;
                if (root2 && p == root) continue;
                const int e = static_cast<int>(es.size());
                es.push_back({F.supp[off + v], F.blen[off + v]});
                ua[v].push_back({p, e});
                ua[p].push_back({v, e});
            }
            if (root2) {
                const int c0 = F.children[F.child_off[off + root]] - off, c1 = F.children[F.child_off[off + root] + 1] - off;
                const int e = static_cast<int>(es.size());
                es.push_back({combine_support(F.supp[off + c0], F.supp[off + c1]), combine_length(F.blen[off + c0], F.blen[off + c1])});
                ua[c0].push_back({c1, e});
                ua[c1].push_back({c0, e});
            }
            if (es.size() >= 64) throw std::logic_error("brute: too many edges for a bitmask");
            std::vector<uint64_t> pm(static_cast<size_t>(F.n) * F.n, 0);
            for (int v = 0; v < m; ++v) {
                if (F.taxon[off + v] < 0) continue;
                std::vector<std::pair<int, int>> up(m, {-2, -1});   // BFS from this leaf: (parent, edge)
                std::queue<int> q;
                up[v] = {-1, -1};
                q.push(v);
                while (!q.empty()) {
                    const int x = q.front(); q.pop();
                    for (auto [y, e] : ua[x]) if (up[y].first == -2) { up[y] = {x, e}; q.push(y); }
                }
                for (int u = 0; u < m; ++u) {
                    if (F.taxon[off + u] < 0) continue;
                    uint64_t mask = 0;
                    for (int w = u; up[w].first >= 0; w = up[w].first) mask |= uint64_t(1) << up[w].second;
                    pm[static_cast<size_t>(F.taxon[off + v]) * F.n + F.taxon[off + u]] = mask;
                }
            }
            edges.push_back(std::move(es));
            pmask.push_back(std::move(pm));
        }
    }

    uint64_t path(int g, int p, int r) const { return pmask[g][static_cast<size_t>(p) * F.n + r]; }

    // weight of the quartet that tree g displays as pq|rs: the anchors are the nodes where p,q and
    // r,s meet; the central path joins them and the pendant paths lead from the leaves to the anchors
    double quartet_weight(int g, int p, int q, int r, int s) const {
        if (mode == 0) return 1.0;
        const uint64_t Epr = path(g, p, r), Eqs = path(g, q, s), Eps = path(g, p, s), Eqr = path(g, q, r);
        const uint64_t Epq = path(g, p, q), Ers = path(g, r, s);
        const uint64_t central = Epr & Eqs;
        if (central == 0 || central != (Eps & Eqr) || central != (Epr & Eps & Eqr & Eqs)) throw std::logic_error("brute: central path");
        const uint64_t pend = (Epq & Epr) | (Epq & Eqr) | (Ers & Epr) | (Ers & Eps);
        if ((pend & central) != 0 || ((Epq & Epr) | central | (Ers & Epr)) != Epr) throw std::logic_error("brute: pendant paths");
        double w = 1.0;
        if (mode & 1) {
            double prod = 1.0;
            for (size_t e = 0; e < edges[g].size(); ++e)
                if (central >> e & 1) prod *= 1.0 - (std::isnan(edges[g][e].s) ? 1.0 : edges[g][e].s);
            w *= 1.0 - prod;
        }
        if (mode & 2) {
            double L = 0.0;
            for (size_t e = 0; e < edges[g].size(); ++e)
                if (pend >> e & 1) L += std::isnan(edges[g][e].l) ? 0.0 : edges[g][e].l;
            w *= std::exp(-L);
        }
        if (!(w >= 0.0 && w <= 1.0)) throw std::logic_error("brute: weight out of range");
        return w;
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
                            const int t = (s1 < s2 && s1 < s3) ? 1 : (s2 < s1 && s2 < s3) ? 2 : (s3 < s1 && s3 < s2) ? 3 : 0;
                            const int topo = t == 1 ? 1 : t == 0 ? 0 : -1;
                            const double w = t == 1 ? quartet_weight(g, a, b, c, e) : t == 2 ? quartet_weight(g, a, c, b, e)
                                           : t == 3 ? quartet_weight(g, a, e, b, c) : 1.0;
                            total += weight[a] * weight[b] * weight[c] * weight[e] * (topo + 1 - g_kappa) * w;
                        }
                }
        }
        return total;
    }

    // topology of quartet a,b,c,e in tree g: 1 ab|ce, 2 ac|be, 3 ae|bc, 0 unresolved
    int quartet(int g, int a, int b, int c, int e) const {
        auto D = [&](int x, int y) { return dist[g][static_cast<size_t>(x) * F.n + y]; };
        const int s1 = D(a, b) + D(c, e), s2 = D(a, c) + D(b, e), s3 = D(a, e) + D(b, c);
        if (s1 < s2 && s1 < s3) return 1;
        if (s2 < s1 && s2 < s3) return 2;
        if (s3 < s1 && s3 < s2) return 3;
        return 0;
    }
};

// resolve_polytomies(F): every internal node has two children, each tree keeps its taxa, every quartet
// is resolved, and quartets resolved in F keep their topology. Returns the number of failures.
long check_resolve(const Forest& F, const Brute& brute, long& quartets) {
    const Forest R = resolve_polytomies(F);
    const Brute rb(R);
    long bad = 0;
    for (int g = 0; g < F.k; ++g) {
        for (int v = R.tree_off[g]; v < R.tree_off[g + 1]; ++v) {
            const int deg = R.child_off[v + 1] - R.child_off[v];
            if (deg != 0 && deg != 2) { std::printf("RESOLVE tree %d: node with %d children\n", g, deg); ++bad; }
        }
        std::vector<int> P;
        for (int t = 0; t < F.n; ++t) {
            if ((F.leaf(g, t) >= 0) != (R.leaf(g, t) >= 0)) { std::printf("RESOLVE tree %d: taxon %d lost or added\n", g, t); ++bad; }
            if (F.leaf(g, t) >= 0) P.push_back(t);
        }
        const int p = static_cast<int>(P.size());
        for (int i = 0; i < p; ++i)
            for (int j = i + 1; j < p; ++j)
                for (int k = j + 1; k < p; ++k)
                    for (int l = k + 1; l < p; ++l) {
                        const int before = brute.quartet(g, P[i], P[j], P[k], P[l]), after = rb.quartet(g, P[i], P[j], P[k], P[l]);
                        if (after == 0 || (before != 0 && before != after)) {
                            std::printf("RESOLVE tree %d: quartet %d %d %d %d: %d -> %d\n", g, P[i], P[j], P[k], P[l], before, after);
                            ++bad;
                        }
                        ++quartets;
                    }
    }
    return bad;
}

double rel(double x, double y) { return std::fabs(x - y) / std::max(1.0, std::fabs(y)); }

// The parser and make_forest keep supports and lengths (a chain through a single-child node is one
// branch), and the reference weights a degree-2 root's two branches as one edge.
void check_attributes() {
    auto eq = [](double x, double y) { return std::fabs(x - y) < 1e-12; };
    Taxa taxa;
    for (const char* l : {"a", "b", "c", "d", "e"}) taxa.add(l);
    const Forest F = make_forest({parse_newick("((a:1,b:2)0.5:3,((c:4,d:5)0.25:6)0.75:1.5,e:7);")}, taxa);
    const int a = F.leaf(0, 0), b = F.leaf(0, 1), c = F.leaf(0, 2), e = F.leaf(0, 4);
    const int ab = F.parent[a], cd = F.parent[c], root = F.parent[ab];
    bool ok = F.nodes() == 8 && F.parent[cd] == root && F.parent[e] == root && F.parent[root] == -1
        && eq(F.blen[a], 1) && eq(F.blen[b], 2) && eq(F.blen[c], 4) && eq(F.blen[e], 7) && std::isnan(F.supp[a])
        && eq(F.supp[ab], 0.5) && eq(F.blen[ab], 3) && eq(F.supp[cd], 0.75) && eq(F.blen[cd], 7.5)
        && std::isnan(F.supp[root]) && std::isnan(F.blen[root]);
    if (!ok) { std::printf("ATTRIBUTES: parsed supports or lengths are wrong\n"); std::exit(1); }
    Taxa t4;
    for (const char* l : {"a", "b", "c", "d"}) t4.add(l);
    const Forest G = make_forest({parse_newick("((a:1,b:2)0.5:3,(c:4,d:5)0.25:6);")}, t4);
    const Brute bs(G, 1), bl(G, 2), bh(G, 3);
    const std::vector<double> ones(4, 1.0);
    const std::vector<int> none(4, -1);
    const std::vector<int> sat{0, 0, 1, 1}, viol{0, 1, 0, 1};   // ab|cd: central edge support 0.5, pendant length 12
    ok = eq(bs.score(sat, ones, none), 0.5) && eq(bl.score(sat, ones, none), std::exp(-12.0))
        && eq(bh.score(sat, ones, none), 0.5 * std::exp(-12.0)) && eq(bs.score(viol, ones, none), -0.5);
    if (!ok) { std::printf("ATTRIBUTES: brute-force weights are wrong\n"); std::exit(1); }
}

}  // namespace

int main(int argc, char** argv) {
    const uint64_t seed = argc > 1 ? std::stoull(argv[1]) : 1;
    const int trials = argc > 2 ? std::stoi(argv[2]) : 300;
    g_kappa = argc > 3 ? std::stod(argv[3]) : 1.0;
    g_gamma = argc > 4 ? std::stod(argv[4]) : 0.0;
    g_weight = argc > 5 ? std::stoi(argv[5]) : 0;
    check_attributes();
    rng.seed(seed);
    double worst_score = 0, worst_gain = 0;
    long nodes_compressed = 0, nodes_plain = 0, resolve_bad = 0, resolve_quartets = 0, weight_bad = 0, weight_checks = 0;
    int checks = 0, child_checks = 0;

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
        const Brute brute(F, g_weight);
        resolve_bad += check_resolve(F, brute, resolve_quartets);
        if (g_weight) {
            // the reference's weights reduce to 1 when every support and length is absent, and are
            // consistent (quartet_weight throws on a broken path identity or a weight outside [0, 1])
            Forest Fn = F;
            std::fill(Fn.supp.begin(), Fn.supp.end(), std::nan(""));
            std::fill(Fn.blen.begin(), Fn.blen.end(), std::nan(""));
            const Brute neutral(Fn, g_weight), plain(F, 0);
            std::vector<int> side(n);
            for (auto& x : side) x = rnd(2);
            const std::vector<double> ones(n, 1.0);
            const std::vector<int> none(n, -1);
            if (neutral.score(side, ones, none) != plain.score(side, ones, none)) {
                std::printf("WEIGHT trial %d: neutral attributes do not give the plain score\n", trial);
                ++weight_bad;
            }
            brute.score(side, ones, none);         // runs the path identities on every resolved quartet
            ++weight_checks;
        }
        auto make = [&](const Forest& f, bool compress) -> std::unique_ptr<IScorer> {
            if (g_weight == 0) return std::make_unique<Scorer>(f, compress);
            WScorer::Options o;
            o.support = g_weight & 1;
            o.length = g_weight & 2;
            return std::make_unique<WScorer>(f, compress, o);
        };

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
            for (int r : d.reals) { weight[r] = dummy_weight(c); dummy_of[r] = top; }
            for (int s : d.subs) self(self, s, c, top);
        };
        for (size_t j = 0; j < sp.dummies.size(); ++j) assign(assign, sp.dummies[j], 1.0, static_cast<int>(j));
        auto sides = [&](const std::vector<int8_t>& pr, const std::vector<int8_t>& pd) {
            std::vector<int> side(n);
            for (size_t i = 0; i < sp.reals.size(); ++i) side[sp.reals[i]] = pr[i];
            for (int t = 0; t < n; ++t) if (dummy_of[t] >= 0) side[t] = pd[dummy_of[t]];
            return side;
        };

        for (bool compress : {true, false}) {
            std::unique_ptr<IScorer> scp = make(F, compress);
            IScorer& sc = *scp;
            sc.set_kappa(g_kappa); sc.set_dummy_gamma(g_gamma);
            sc.setup(sp, pool);
            (compress ? nodes_compressed : nodes_plain) += sc.forest().nodes();
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
                if (e > 1e-9) std::printf("REAL GAIN trial %d unit %zu (compress %d): %.10g vs brute %.10g\n", trial, i, compress, gr[i], expect);
            }
            for (size_t j = 0; j < sp.dummies.size(); ++j) {
                auto pd = sp.part_d; pd[j] ^= 1;
                const double expect = brute.score(sides(sp.part_r, pd), weight, dummy_of) - base;
                const double e = rel(gd[j], expect);
                worst_gain = std::max(worst_gain, e);
                if (e > 1e-9) std::printf("DUMMY GAIN trial %d dummy %zu (compress %d): %.10g vs brute %.10g\n", trial, j, compress, gd[j], expect);
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
            // a batch of bipartitions (resets and real-taxon moves): against brute force, and bit for bit
            // against the same steps through set_partition / swap_real / score
            ScoreBatch batch;
            std::vector<double> step;
            std::unique_ptr<IScorer> seqp = make(F, compress);
            IScorer& seq = *seqp;
            seq.set_kappa(g_kappa); seq.set_dummy_gamma(g_gamma);
            seq.setup(sp, pool);
            for (int q = 0; q < 6; ++q) {
                const bool reset = q == 0 || rnd(3) == 0;
                batch.reset.push_back(reset);
                if (reset) {
                    for (auto& x : pr) x = static_cast<int8_t>(rnd(2));
                    for (auto& x : pd) x = static_cast<int8_t>(rnd(2));
                    seq.set_partition(pr, pd);
                } else {
                    for (int i = 0; i < static_cast<int>(pr.size()); ++i)
                        if (rnd(4) == 0) { pr[i] ^= 1; seq.swap_real(i); batch.moved.push_back(i); }
                }
                batch.off.push_back(static_cast<int>(batch.moved.size()));
                batch.sides.insert(batch.sides.end(), pr.begin(), pr.end());
                batch.sides.insert(batch.sides.end(), pd.begin(), pd.end());
                step.push_back(seq.score());
                worst_score = std::max(worst_score, rel(step.back(), brute.score(sides(pr, pd), weight, dummy_of)));
            }
            std::vector<double> out;
            sc.scores(batch, out);
            for (int q = 0; q < batch.size(); ++q)
                if (out[q] != step[q]) {
                    std::printf("BATCH trial %d item %d (compress %d): %.17g vs step by step %.17g\n", trial, q, compress, out[q], step[q]);
                    worst_score = std::max(worst_score, 1.0);
                }
            // afterwards the scorer holds the last bipartition
            std::vector<double> gr2, gd2;
            sc.gains(gr, gd);
            seq.gains(gr2, gd2);
            bool same = gr == gr2 && gd == gd2;
            for (size_t i = 0; i < pr.size(); ++i) same = same && sc.part_real(static_cast<int>(i)) == pr[i];
            for (size_t j = 0; j < pd.size(); ++j) same = same && sc.part_dummy(static_cast<int>(j)) == pd[j];
            if (!same) {
                std::printf("BATCH trial %d (compress %d): state after the batch differs\n", trial, compress);
                worst_gain = std::max(worst_gain, 1.0);
            }
        }
        // child subproblems built from the parent's forest, as in the divide step (Scorer::setup_child)
        {
            std::unique_ptr<IScorer> parentp = make(F, true);
            IScorer& parent = *parentp;
            parent.set_kappa(g_kappa); parent.set_dummy_gamma(g_gamma);
            parent.setup(sp, pool);
            parent.set_partition(sp.part_r, sp.part_d);
            for (int side = 0; side < 2; ++side) {
                Subproblem ch;
                std::vector<int> other_r, other_d;
                for (size_t i = 0; i < sp.reals.size(); ++i) (sp.part_r[i] == side ? ch.reals : other_r).push_back(sp.reals[i]);
                for (size_t j = 0; j < sp.dummies.size(); ++j) (sp.part_d[j] == side ? ch.dummies : other_d).push_back(sp.dummies[j]);
                if (other_r.empty() && other_d.empty()) continue;
                ch.dummies.push_back(pool.add(other_r, other_d));
                for (size_t i = 0; i < ch.reals.size(); ++i) ch.part_r.push_back(static_cast<int8_t>(rnd(2)));
                for (size_t j = 0; j < ch.dummies.size(); ++j) ch.part_d.push_back(static_cast<int8_t>(rnd(2)));
                std::vector<double> cw(n, 1.0);
                std::vector<int> cdum(n, -1);
                auto cassign = [&](auto&& self, int id, double denom, int top) -> void {
                    const Dummy& d = pool[id];
                    const double c = denom * static_cast<double>(d.reals.size() + d.subs.size());
                    for (int r : d.reals) { cw[r] = dummy_weight(c); cdum[r] = top; }
                    for (int x : d.subs) self(self, x, c, top);
                };
                for (size_t j = 0; j < ch.dummies.size(); ++j) cassign(cassign, ch.dummies[j], 1.0, static_cast<int>(j));
                auto csides = [&](const std::vector<int8_t>& pr, const std::vector<int8_t>& pd) {
                    std::vector<int> sd(n);
                    for (size_t i = 0; i < ch.reals.size(); ++i) sd[ch.reals[i]] = pr[i];
                    for (int t = 0; t < n; ++t) if (cdum[t] >= 0) sd[t] = pd[cdum[t]];
                    return sd;
                };
                std::unique_ptr<IScorer> derivedp = make(F, true);
                IScorer& derived = *derivedp;
                derived.set_kappa(g_kappa); derived.set_dummy_gamma(g_gamma);
                derived.setup_child(ch, parent, side);
                derived.set_partition(ch.part_r, ch.part_d);
                std::vector<double> cgr, cgd;
                const double fast = derived.gains(cgr, cgd);
                const double base = brute.score(csides(ch.part_r, ch.part_d), cw, cdum);
                worst_score = std::max(worst_score, rel(fast, base));
                for (size_t i = 0; i < ch.reals.size(); ++i) {
                    auto pr = ch.part_r; pr[i] ^= 1;
                    const double e = rel(cgr[i], brute.score(csides(pr, ch.part_d), cw, cdum) - base);
                    worst_gain = std::max(worst_gain, e);
                    if (e > 1e-9) std::printf("CHILD REAL GAIN trial %d side %d unit %zu\n", trial, side, i);
                }
                for (size_t j = 0; j < ch.dummies.size(); ++j) {
                    auto pd = ch.part_d; pd[j] ^= 1;
                    const double e = rel(cgd[j], brute.score(csides(ch.part_r, pd), cw, cdum) - base);
                    worst_gain = std::max(worst_gain, e);
                    if (e > 1e-9) std::printf("CHILD DUMMY GAIN trial %d side %d dummy %zu\n", trial, side, j);
                }
                ++child_checks;
            }
        }
        ++checks;
    }
    std::printf("brute_check seed=%llu subproblems=%d (+%d derived children)  max rel err: score %.2e, gains %.2e  (forest nodes: compressed %ld, plain %ld)\n",
                static_cast<unsigned long long>(seed), checks, child_checks, worst_score, worst_gain, nodes_compressed, nodes_plain);
    std::printf("resolve_polytomies: %ld quartets checked, %ld failures\n", resolve_quartets, resolve_bad);
    if (g_weight) std::printf("weighted reference (mode %d): %ld forests checked, %ld failures\n", g_weight, weight_checks, weight_bad);

    // quartet_support vs enumeration: random species trees (some with polytomies) against random gene
    // trees with polytomies and missing taxa
    long qs_bad = 0, qs_quartets = 0;
    for (int trial = 0; trial < trials; ++trial) {
        const int n = 4 + rnd(12), k = 1 + rnd(5);
        Taxa taxa;
        std::vector<std::string> labels;
        for (int i = 0; i < n; ++i) { labels.push_back("t" + std::to_string(i)); taxa.add(labels.back()); }
        std::vector<ParsedTree> trees;
        for (int g = 0; g < k; ++g) {
            std::vector<std::string> present = labels;
            if (unif() < 0.5) {
                const int drop = rnd(std::max(1, n - 3));
                for (int i = 0; i < drop; ++i) present.erase(present.begin() + rnd(static_cast<int>(present.size())));
            }
            trees.push_back(parse_newick(random_tree(present, unif() < 0.5 ? 0.0 : 0.4, unif() < 0.3)));
        }
        const Forest F = make_forest(trees, taxa);
        const Forest S = make_forest({parse_newick(random_tree(labels, unif() < 0.7 ? 0.0 : 0.4, unif() < 0.5))}, taxa);
        const Brute bg(F), bs(S);
        QuartetCounts want;
        for (int g = 0; g < F.k; ++g) {
            std::vector<int> P;
            for (int t = 0; t < n; ++t) if (F.leaf(g, t) >= 0) P.push_back(t);
            const int p = static_cast<int>(P.size());
            for (int i = 0; i < p; ++i)
                for (int j = i + 1; j < p; ++j)
                    for (int a = j + 1; a < p; ++a)
                        for (int b = a + 1; b < p; ++b) {
                            const int tg = bg.quartet(g, P[i], P[j], P[a], P[b]);
                            if (tg == 0) continue;
                            ++want.resolved;
                            want.agree += tg == bs.quartet(0, P[i], P[j], P[a], P[b]);
                        }
        }
        const QuartetCounts got = quartet_support(F, S, 1 + trial % 3);
        if (got.agree != want.agree || got.resolved != want.resolved) {
            std::printf("QSCORE trial %d: agree %lld (want %lld), resolved %lld (want %lld)\n", trial,
                        static_cast<long long>(got.agree), static_cast<long long>(want.agree),
                        static_cast<long long>(got.resolved), static_cast<long long>(want.resolved));
            ++qs_bad;
        }
        qs_quartets += want.resolved;
    }
    std::printf("quartet_support: %d species trees, %ld resolved gene quartets, %ld failures\n", trials, qs_quartets, qs_bad);
    return (worst_score > 1e-9 || worst_gain > 1e-9 || resolve_bad > 0 || qs_bad > 0 || weight_bad > 0) ? 1 : 0;
}
