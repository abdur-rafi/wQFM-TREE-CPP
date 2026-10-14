#include "consensus.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "trace.hpp"

namespace wqfm {
namespace {

void erase_first(std::vector<int>& v, int x) {
    auto it = std::find(v.begin(), v.end(), x);
    if (it != v.end()) v.erase(it);
}

}  // namespace

ConsensusPartitioner::ConsensusPartitioner(const ParsedTree& t, const Taxa& taxa)
    : childs_(t.children), parent_(t.parent), root_(t.root), n_(taxa.size()) {
    taxon_.assign(t.parent.size(), -1);
    std::vector<char> seen(n_, 0);
    for (size_t v = 0; v < t.parent.size(); ++v) {
        if (!t.children[v].empty()) continue;
        const int id = taxa.find(t.label[v]);
        if (id < 0) throw std::runtime_error("consensus tree: unknown taxon '" + t.label[v] + "'");
        if (seen[id]) throw std::runtime_error("consensus tree: taxon '" + t.label[v] + "' appears twice");
        seen[id] = 1;
        taxon_[v] = id;
    }
    for (int id = 0; id < n_; ++id)
        if (!seen[id]) throw std::runtime_error("consensus tree is missing taxon '" + taxa.label[id] + "'");

    if (childs_[root_].size() > 2) balance_root();

    // post-order in child-list order (Java Tree.topSort)
    std::vector<std::pair<int, size_t>> st{{root_, 0}};
    while (!st.empty()) {
        auto& [v, k] = st.back();
        if (k < childs_[v].size()) {
            const int c = childs_[v][k++];
            st.emplace_back(c, 0);
        } else {
            topo_.push_back(v);
            st.pop_back();
        }
    }
}

// Java Tree.balanceRoot: re-root next to the node whose subtree node count is closest to half
void ConsensusPartitioner::balance_root() {
    const int n = static_cast<int>(parent_.size());
    std::vector<int> cnt(n, 1);
    for (int v = 0; v < n; ++v)                   // creation order: children before parents
        for (int c : childs_[v]) cnt[v] += cnt[c];
    int closest = root_, diff = n;
    for (int i = 0; i < n; ++i) {
        const int v = std::abs(n / 2 - cnt[i]);
        if (v < diff) { diff = v; closest = i; }
    }
    const int cp = parent_[closest];
    if (cp < 0) throw std::runtime_error("consensus tree: cannot re-balance the root");
    erase_first(childs_[cp], closest);
    int curr = cp, currp = parent_[curr];
    while (curr != -1 && currp != -1) {
        childs_[curr].push_back(currp);
        erase_first(childs_[currp], curr);
        const int temp = currp;
        currp = parent_[currp];
        parent_[temp] = curr;
        curr = temp;
    }
    const int r = n;
    childs_.push_back({closest, cp});
    parent_.push_back(-1);
    taxon_.push_back(-1);
    parent_[closest] = r;
    parent_[cp] = r;
    root_ = r;
}

void ConsensusPartitioner::make(Subproblem& sp, bool allow_singleton, Scorer& scorer, const DummyPool& pool,
                                Trace* trace, const Replay* replay) const {
    const int R = static_cast<int>(sp.reals.size()), D = static_cast<int>(sp.dummies.size());
    const std::string key = (trace || replay) ? sub_key(sp, pool) : std::string();   // O(n): only when needed
    if (trace) trace->line("P " + key);

    // weights as in Java: 1 for real taxa, Eq. 1 for taxa under a dummy
    std::vector<double> weight(n_, 0.0);
    std::vector<int> real_index(n_, -1), in_dummy(n_, -1);
    for (int i = 0; i < R; ++i) { weight[sp.reals[i]] = 1; real_index[sp.reals[i]] = i; }
    for (int j = 0; j < D; ++j) {
        pool.coeffs(sp.dummies[j], 1.0, weight);
        for (int y : pool[sp.dummies[j]].flat) { weight[y] = 1 / weight[y]; in_dummy[y] = j; }
    }

    const int N = static_cast<int>(childs_.size());
    std::vector<int> real_cnt(N, 0);
    std::vector<double> dw(static_cast<size_t>(N) * D, 0.0);

    std::vector<int8_t> pr(R), pd(D), best_r, best_d;
    std::unordered_map<std::string, int> cache;   // bipartition -> its index in the batch
    std::vector<int> stack;
    std::vector<std::pair<int, int>> cands;  // (subtree node, batch index) in Java's order
    Scorer::Batch batch;                     // the distinct bipartitions, scored together below
    std::vector<int8_t> cur_r, cur_d;        // the batch's last bipartition
    bool fresh = true;                       // batch empty: the scorer holds no partition yet
    const SubForest& sf = scorer.forest();
    const int max_real_moves = std::max(1, sf.k ? sf.nodes() / sf.k / 8 : 1);   // ~ cost of one rebuild

    // the bipartition of candidate node c: its subtree's units on side 1 (a dummy where it has >= half
    // its weight)
    auto candidate = [&](int c, std::vector<int8_t>& part_r, std::vector<int8_t>& part_d) {
        part_r.assign(R, 0);
        stack.assign(1, c);
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            if (childs_[u].empty()) { if (real_index[taxon_[u]] >= 0) part_r[real_index[taxon_[u]]] = 1; }
            else stack.insert(stack.end(), childs_[u].begin(), childs_[u].end());
        }
        part_d.resize(D);
        const double* dc = &dw[static_cast<size_t>(c) * D];
        for (int j = 0; j < D; ++j) part_d[j] = dc[j] >= .5 ? 1 : 0;
    };

    for (int v : topo_) {
        double* dv = &dw[static_cast<size_t>(v) * D];
        if (childs_[v].empty()) {
            const int t = taxon_[v];
            if (real_index[t] >= 0) real_cnt[v] = 1;
            else if (weight[t] != 0) dv[in_dummy[t]] = weight[t];
            continue;
        }
        for (int c : childs_[v]) {
            const double* dc = &dw[static_cast<size_t>(c) * D];
            int pa = real_cnt[c], pb = R - pa;
            for (int j = 0; j < D; ++j) {
                dv[j] += dc[j];
                if (dc[j] >= .5) ++pa; else ++pb;
            }
            real_cnt[v] += real_cnt[c];
            if (!((pa > 1 && pb > 1) || (allow_singleton && pa >= 1 && pb >= 1))) continue;

            candidate(c, pr, pd);

            const std::string bits = part_bits(pr, pd);
            auto it = cache.find(bits);
            if (it == cache.end()) {
                // reach this candidate from the previous one: consecutive candidates usually differ in a
                // few real taxa, and moving those is cheaper than rebuilding (a dummy move costs a full
                // pass anyway)
                int diff_r = 0, diff_d = 0;
                if (!fresh) {
                    for (int i = 0; i < R; ++i) diff_r += pr[i] != cur_r[i];
                    for (int j = 0; j < D; ++j) diff_d += pd[j] != cur_d[j];
                }
                const bool reset = fresh || diff_d > 0 || diff_r > max_real_moves;
                batch.reset.push_back(reset);
                if (!reset)
                    for (int i = 0; i < R; ++i) if (pr[i] != cur_r[i]) batch.moved.push_back(i);
                batch.off.push_back(static_cast<int>(batch.moved.size()));
                batch.sides.insert(batch.sides.end(), pr.begin(), pr.end());
                batch.sides.insert(batch.sides.end(), pd.begin(), pd.end());
                cur_r = pr;
                cur_d = pd;
                fresh = false;
                it = cache.emplace(bits, batch.size() - 1).first;
            }
            cands.emplace_back(c, it->second);
        }
    }
    std::vector<double> score;
    scorer.scores(batch, score);
    if (trace) for (const auto& [c, q] : cands) trace->line("C " + Trace::num(score[q]));

    // Java keeps the first candidate with the largest score (strict >). Its tied scores are bit-identical
    // (the same quartet weights summed), ours carry rounding noise that would pick an arbitrary one of
    // them, so scores within a tolerance of the largest count as tied and the first of them wins. Noise
    // is below ~1e-13 of the largest |score| (or of the weight scale when all scores are ~0); distinct
    // scores of the best candidates differ by more than 1e-4 of it on the ASTRAL-II and simulated data.
    if (!cands.empty()) {
        double max_s = score[cands[0].second], mag = 0;
        for (const auto& [c, q] : cands) {
            max_s = std::max(max_s, score[q]);
            mag = std::max(mag, std::fabs(score[q]));
        }
        const double tol = 1e-9 * mag + 1e-16 * scorer.scale();
        for (const auto& [c, q] : cands) {
            if (score[q] >= max_s - tol) {
                candidate(c, best_r, best_d);
                break;
            }
        }
    } else {                // Java falls back to RandPartition: first half of the units on side 0
        int pa = (R + D) / 2;
        best_r.assign(R, 1);
        best_d.assign(D, 1);
        for (int i = 0; i < R && pa > 0; ++i, --pa) best_r[i] = 0;
        for (int j = 0; j < D && pa > 0; ++j, --pa) best_d[j] = 0;
    }
    if (replay) {
        auto it = replay->parts.find(key);
        if (it != replay->parts.end() && !it->second.result.empty()) {
            auto [r, d] = parse_bits(it->second.result);
            if (static_cast<int>(r.size()) != R || static_cast<int>(d.size()) != D)
                throw std::runtime_error("replay: partition size mismatch for " + key);
            best_r = std::move(r);
            best_d = std::move(d);
        }
    }
    sp.part_r = best_r;
    sp.part_d = best_d;
    if (trace) trace->line("R " + part_bits(sp.part_r, sp.part_d));
}

}  // namespace wqfm
