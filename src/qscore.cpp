#include "qscore.hpp"

#include <algorithm>
#include <atomic>
#include <exception>
#include <stdexcept>
#include <thread>
#include <vector>

namespace wqfm {
namespace {

inline int64_t pairs(int64_t x) { return x * (x - 1) / 2; }

// The subtrees of a node are those of its children and the rest of the tree. A quartet ab|cd that a
// tree resolves has two join nodes: one where a and b lie in different subtrees and c, d in a third,
// and one where c and d do. With I[i][p] = |species subtree i of u ∩ gene subtree p of w|, the sum
// over species nodes u, gene nodes w, species subtrees k and i < j (neither k) and gene subtrees r of
//     C(I[k][r], 2) * #{a in i, b in j: a and b in different gene subtrees, neither r}
// counts every quartet that both trees resolve the same way twice (once per pair of join nodes).
// The quartets a gene tree resolves are counted twice the same way, from the gene tree alone.
// C[w * mS + u] = |gene clade w ∩ species clade u|.
void count_tree(const Forest& G, int g, const Forest& S, const std::vector<int>& inner, std::vector<int32_t>& C,
                std::vector<int64_t>& I, std::vector<int64_t>& z, QuartetCounts& out) {
    const int off = G.tree_off[g], mG = G.tree_off[g + 1] - off, mS = S.nodes(), rootS = mS - 1;
    C.resize(static_cast<size_t>(mG) * mS);
    for (int w = 0; w < mG; ++w) {                 // children precede their parent
        int32_t* row = &C[static_cast<size_t>(w) * mS];
        std::fill(row, row + mS, 0);
        if (G.taxon[off + w] >= 0) {
            for (int u = S.leaf(0, G.taxon[off + w]); u >= 0; u = S.parent[u]) row[u] = 1;
            continue;
        }
        for (int i = G.child_off[off + w]; i < G.child_off[off + w + 1]; ++i) {
            const int32_t* cr = &C[static_cast<size_t>(G.children[i] - off) * mS];
            for (int u = 0; u < mS; ++u) row[u] += cr[u];
        }
    }
    const int32_t* top = &C[static_cast<size_t>(mG - 1) * mS];   // species clades restricted to the gene's taxa
    const int64_t n = top[rootS];
    int64_t agree2 = 0, resolved2 = 0;
    for (int w = 0; w < mG; ++w) {
        const int cb = G.child_off[off + w], ce = G.child_off[off + w + 1];
        if (cb == ce) continue;
        const int P = ce - cb + 1;                 // gene subtrees: the children, then the rest
        const int32_t* row_w = &C[static_cast<size_t>(w) * mS];
        z.assign(P, 0);
        int64_t T = 0, Q = 0;
        for (int p = 0; p + 1 < P; ++p) z[p] = C[static_cast<size_t>(G.children[cb + p] - off) * mS + rootS];
        z[P - 1] = n - row_w[rootS];
        for (int p = 0; p < P; ++p) { T += z[p]; Q += z[p] * z[p]; }
        for (int r = 0; r < P; ++r) resolved2 += pairs(z[r]) * (((T - z[r]) * (T - z[r]) - (Q - z[r] * z[r])) / 2);

        for (int u : inner) {
            const int sb = S.child_off[u], se = S.child_off[u + 1], E = se - sb + 1;   // species subtrees
            I.resize(static_cast<size_t>(E) * P);
            for (int i = 0; i + 1 < E; ++i) {
                const int s = S.children[sb + i];
                for (int p = 0; p + 1 < P; ++p) I[i * P + p] = C[static_cast<size_t>(G.children[cb + p] - off) * mS + s];
                I[i * P + P - 1] = top[s] - row_w[s];
            }
            for (int p = 0; p + 1 < P; ++p) {
                const int32_t* rc = &C[static_cast<size_t>(G.children[cb + p] - off) * mS];
                I[(E - 1) * P + p] = rc[rootS] - rc[u];
            }
            I[(E - 1) * P + P - 1] = n - top[u] - row_w[rootS] + row_w[u];

            for (int i = 0; i < E; ++i)
                for (int j = i + 1; j < E; ++j) {
                    const int64_t* Ii = &I[i * P];
                    const int64_t* Ij = &I[j * P];
                    int64_t Ai = 0, Aj = 0, Pij = 0;
                    for (int p = 0; p < P; ++p) { Ai += Ii[p]; Aj += Ij[p]; Pij += Ii[p] * Ij[p]; }
                    if (Ai == 0 || Aj == 0) continue;
                    for (int k = 0; k < E; ++k) {
                        if (k == i || k == j) continue;
                        const int64_t* Ik = &I[k * P];
                        for (int r = 0; r < P; ++r) {
                            if (Ik[r] < 2) continue;
                            agree2 += pairs(Ik[r]) * ((Ai - Ii[r]) * (Aj - Ij[r]) - (Pij - Ii[r] * Ij[r]));
                        }
                    }
                }
        }
    }
    out.agree += agree2 / 2;
    out.resolved += resolved2 / 2;
}

}  // namespace

QuartetCounts quartet_support(const Forest& genes, const Forest& species, int threads) {
    if (species.k != 1) throw std::runtime_error("quartet_support: the species forest must hold one tree");
    for (int v = 0; v < genes.nodes(); ++v)
        if (genes.taxon[v] >= 0 && species.leaf(0, genes.taxon[v]) < 0)
            throw std::runtime_error("the species tree is missing a taxon of the gene trees");
    std::vector<int> inner;
    for (int u = 0; u < species.nodes(); ++u)
        if (species.child_off[u + 1] > species.child_off[u]) inner.push_back(u);

    threads = std::max(1, std::min(threads, genes.k));
    std::vector<QuartetCounts> part(threads);
    std::vector<std::exception_ptr> err(threads);
    std::atomic<int> next{0};
    auto work = [&](int id) {
        try {
            std::vector<int32_t> C;
            std::vector<int64_t> I, z;
            for (int g; (g = next.fetch_add(1)) < genes.k;) count_tree(genes, g, species, inner, C, I, z, part[id]);
        } catch (...) {
            err[id] = std::current_exception();
        }
    };
    std::vector<std::thread> pool;
    for (int id = 1; id < threads; ++id) pool.emplace_back(work, id);
    work(0);
    for (auto& t : pool) t.join();
    for (auto& e : err) if (e) std::rethrow_exception(e);
    QuartetCounts total;
    for (const auto& p : part) { total.agree += p.agree; total.resolved += p.resolved; }
    return total;
}

}  // namespace wqfm
