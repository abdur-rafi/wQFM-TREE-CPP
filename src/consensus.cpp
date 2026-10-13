#include "consensus.hpp"

#include <algorithm>
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
    const std::string key = sub_key(sp, pool);
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
    scorer.setup(sp, pool);

    std::vector<int8_t> pr(R), pd(D), best_r, best_d;
    std::unordered_map<std::string, double> cache;
    std::vector<int> stack;
    int min_node = -1;
    double max_score = 0;

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

            std::fill(pr.begin(), pr.end(), 0);
            stack.assign(1, c);
            while (!stack.empty()) {
                const int u = stack.back();
                stack.pop_back();
                if (childs_[u].empty()) { if (real_index[taxon_[u]] >= 0) pr[real_index[taxon_[u]]] = 1; }
                else stack.insert(stack.end(), childs_[u].begin(), childs_[u].end());
            }
            for (int j = 0; j < D; ++j) pd[j] = dc[j] >= .5 ? 1 : 0;

            const std::string bits = part_bits(pr, pd);
            auto it = cache.find(bits);
            double s;
            if (it != cache.end()) {
                s = it->second;
            } else {
                scorer.set_partition(pr, pd);
                s = scorer.score();
                cache.emplace(bits, s);
            }
            if (trace) trace->line("C " + Trace::num(s));
            if (min_node == -1 || s > max_score) {
                max_score = s;
                min_node = c;
                best_r = pr;
                best_d = pd;
            }
        }
    }

    if (min_node == -1) {   // Java falls back to RandPartition: first half of the units on side 0
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
