#include "forest.hpp"

#include <stdexcept>

namespace wqfm {

int Taxa::add(const std::string& l) {
    auto [it, inserted] = index.emplace(l, static_cast<int>(label.size()));
    if (inserted) label.push_back(l);
    return it->second;
}

int Taxa::find(const std::string& l) const {
    auto it = index.find(l);
    return it == index.end() ? -1 : it->second;
}

Forest make_forest(const std::vector<ParsedTree>& trees, const Taxa& taxa) {
    Forest f;
    f.n = taxa.size();
    f.k = static_cast<int>(trees.size());
    f.tree_off.push_back(0);
    f.child_off.push_back(0);
    f.leaf_of.assign(static_cast<size_t>(f.k) * f.n, -1);
    std::vector<int> rep, idx;
    for (int g = 0; g < f.k; ++g) {
        const ParsedTree& t = trees[g];
        const int off = f.tree_off.back();
        const int m = static_cast<int>(t.parent.size());
        if (t.root != m - 1) throw std::runtime_error("internal: root is not the last node");
        // contract single-child nodes (e.g. "(a,(b),c)"): same topology, and every internal node then
        // has at least two children (the Java version exits on such nodes)
        rep.assign(m, -1);
        idx.assign(m, -1);
        int kept = 0;
        for (int v = 0; v < m; ++v) {
            if (t.children[v].size() == 1) { rep[v] = rep[t.children[v][0]]; continue; }
            rep[v] = v;
            idx[v] = kept++;
        }
        if (idx[rep[t.root]] != kept - 1) throw std::runtime_error("internal: contracted root is not the last node");
        const int base = static_cast<int>(f.parent.size());
        f.parent.resize(base + kept, -1);
        for (int v = 0; v < m; ++v) {
            if (idx[v] < 0) continue;
            for (int c : t.children[v]) {
                f.children.push_back(off + idx[rep[c]]);
                f.parent[off + idx[rep[c]]] = off + idx[v];
            }
            f.child_off.push_back(static_cast<int>(f.children.size()));
            int tax = -1;
            if (t.children[v].empty()) {
                tax = taxa.find(t.label[v]);
                if (tax < 0) throw std::runtime_error("unknown taxon '" + t.label[v] + "' in gene tree " + std::to_string(g + 1));
                int& slot = f.leaf_of[static_cast<size_t>(g) * f.n + tax];
                if (slot >= 0) throw std::runtime_error("taxon '" + t.label[v] + "' appears twice in gene tree " + std::to_string(g + 1));
                slot = off + idx[v];
            }
            f.taxon.push_back(tax);
        }
        f.tree_off.push_back(off + kept);
    }
    return f;
}

}  // namespace wqfm
