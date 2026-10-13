#include "greedy.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace wqfm {
namespace {

struct Key {
    uint64_t a, b;
    bool operator==(const Key& o) const { return a == o.a && b == o.b; }
};
struct KeyHash {
    size_t operator()(const Key& k) const { return static_cast<size_t>(k.a ^ (k.b * 0x9E3779B97F4A7C15ULL)); }
};
struct SplitInfo {
    int count = 0;
    int tree = -1, node = -1;   // a gene-tree node defining it
    bool flip = false;          // the clade is the complement of that node's subtree
    int order = 0;              // first appearance
};

uint64_t splitmix(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

}  // namespace

std::string greedy_consensus_newick(const Forest& f, const Taxa& taxa) {
    const int n = f.n;
    if (n < 4) throw std::runtime_error("consensus: fewer than 4 taxa");
    for (int g = 0; g < f.k; ++g)
        for (int t = 0; t < n; ++t)
            if (f.leaf(g, t) < 0)
                throw std::runtime_error("the built-in consensus needs complete gene trees (taxon '" + taxa.label[t] +
                                         "' is missing from gene tree " + std::to_string(g + 1) + "); pass a consensus tree with -c");

    // 128-bit XOR hashes of clades; a split is keyed by its side without taxon 0
    std::vector<Key> key(n);
    uint64_t seed = 0x5eed;
    Key all{0, 0};
    for (int t = 0; t < n; ++t) {
        key[t] = {splitmix(seed), splitmix(seed)};
        all.a ^= key[t].a;
        all.b ^= key[t].b;
    }
    std::unordered_map<Key, SplitInfo, KeyHash> splits;
    std::vector<Key> order;
    std::vector<Key> h;
    std::vector<int> size;
    std::vector<char> has0;
    for (int g = 0; g < f.k; ++g) {
        const int off = f.tree_off[g], m = f.tree_off[g + 1] - off, root = off + m - 1;
        h.assign(m, Key{0, 0});
        size.assign(m, 0);
        has0.assign(m, 0);
        for (int v = off; v <= root; ++v) {
            const int lv = v - off;
            if (f.taxon[v] >= 0) {
                h[lv] = key[f.taxon[v]];
                size[lv] = 1;
                has0[lv] = f.taxon[v] == 0;
                continue;
            }
            for (int e = f.child_off[v]; e < f.child_off[v + 1]; ++e) {
                const int c = f.children[e] - off;
                h[lv].a ^= h[c].a;
                h[lv].b ^= h[c].b;
                size[lv] += size[c];
                has0[lv] |= has0[c];
            }
        }
        const bool root2 = f.child_off[root + 1] - f.child_off[root] == 2;
        for (int v = off; v < root; ++v) {
            const int lv = v - off;
            // the two edges below a degree-2 root are one edge of the unrooted tree
            if (root2 && f.parent[v] == root && f.children[f.child_off[root] + 1] == v) continue;
            if (size[lv] < 2 || size[lv] > n - 2) continue;
            const bool flip = has0[lv];
            const Key k = flip ? Key{all.a ^ h[lv].a, all.b ^ h[lv].b} : h[lv];
            auto [it, inserted] = splits.try_emplace(k);
            if (inserted) {
                it->second = SplitInfo{0, g, v, flip, static_cast<int>(order.size())};
                order.push_back(k);
            }
            ++it->second.count;
        }
    }
    std::vector<const SplitInfo*> ranked;
    ranked.reserve(order.size());
    for (const Key& k : order) ranked.push_back(&splits.at(k));
    std::stable_sort(ranked.begin(), ranked.end(), [](const SplitInfo* x, const SplitInfo* y) { return x->count > y->count; });

    // consensus tree rooted where taxon 0 attaches: leaf node t = taxon t (t >= 1), node n = root
    std::vector<int> parent(n + 1, n), csize(n + 1, 1);
    std::vector<std::vector<int>> kids(n + 1);
    parent[n] = -1;
    csize[n] = n - 1;
    for (int t = 1; t < n; ++t) kids[n].push_back(t);

    // nodes: n leaves/root slots plus at most n - 3 accepted splits
    const size_t max_nodes = 2 * static_cast<size_t>(n);
    parent.reserve(max_nodes);
    csize.reserve(max_nodes);
    kids.reserve(max_nodes);
    std::vector<int> clade, stamp(max_nodes, 0), pos(max_nodes, -1), cnt(max_nodes, 0), path, touched, stack;
    std::vector<char> in_clade(n, 0);
    int tag = 0, accepted = 0;
    for (const SplitInfo* sp : ranked) {
        if (accepted == n - 3) break;                      // fully resolved
        // taxa of the split (the side without taxon 0)
        clade.clear();
        stack.assign(1, sp->node);
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            if (f.taxon[v] >= 0) clade.push_back(f.taxon[v]);
            else for (int e = f.child_off[v]; e < f.child_off[v + 1]; ++e) stack.push_back(f.children[e]);
        }
        if (sp->flip) {
            for (int t : clade) in_clade[t] = 1;
            std::vector<int> comp;
            for (int t = 1; t < n; ++t) if (!in_clade[t]) comp.push_back(t);
            for (int t : clade) in_clade[t] = 0;
            clade.swap(comp);
        }
        // lowest common ancestor of the clade in the current consensus tree
        ++tag;
        path.clear();
        for (int v = clade[0]; v != -1; v = parent[v]) { stamp[v] = tag; pos[v] = static_cast<int>(path.size()); path.push_back(v); }
        int top = 0;
        for (size_t i = 1; i < clade.size(); ++i) {
            int v = clade[i];
            while (stamp[v] != tag) v = parent[v];
            top = std::max(top, pos[v]);
        }
        const int lca = path[top];
        if (csize[lca] == static_cast<int>(clade.size())) continue;   // already in the tree
        // every child of the LCA must lie entirely inside or outside the clade
        touched.clear();
        for (int t : clade) {
            int v = t;
            while (parent[v] != lca) v = parent[v];
            if (cnt[v]++ == 0) touched.push_back(v);
        }
        bool ok = true;
        for (int x : touched) ok = ok && cnt[x] == csize[x];
        if (ok) {
            const int u = static_cast<int>(parent.size());
            parent.push_back(lca);
            csize.push_back(static_cast<int>(clade.size()));
            kids.emplace_back();
            std::vector<int>& lk = kids[lca];
            for (int x : touched) {
                lk.erase(std::find(lk.begin(), lk.end(), x));
                parent[x] = u;
                kids[u].push_back(x);
            }
            lk.push_back(u);
            ++accepted;
        }
        for (int x : touched) cnt[x] = 0;
    }

    // Newick, with taxon 0 beside the root's children
    std::string out = "(" + taxa.label[0];
    struct Frame { int v; size_t next; };
    std::vector<Frame> st{{n, 0}};
    while (!st.empty()) {
        Frame& fr = st.back();
        if (fr.next == kids[fr.v].size()) {
            out += ')';
            st.pop_back();
            continue;
        }
        const bool first = fr.next == 0 && fr.v != n;   // the root's children all follow taxon 0
        const int c = kids[fr.v][fr.next++];
        if (!first) out += ',';
        if (kids[c].empty()) out += taxa.label[c];
        else { out += '('; st.push_back({c, 0}); }      // invalidates fr
    }
    return out + ";";
}

}  // namespace wqfm
