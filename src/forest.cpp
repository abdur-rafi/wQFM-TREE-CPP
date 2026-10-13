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
    for (int g = 0; g < f.k; ++g) {
        const ParsedTree& t = trees[g];
        const int off = f.tree_off.back();
        const int m = static_cast<int>(t.parent.size());
        for (int v = 0; v < m; ++v) {
            f.parent.push_back(t.parent[v] < 0 ? -1 : off + t.parent[v]);
            for (int c : t.children[v]) f.children.push_back(off + c);
            f.child_off.push_back(static_cast<int>(f.children.size()));
            int tax = -1;
            if (t.children[v].empty()) {
                tax = taxa.find(t.label[v]);
                if (tax < 0) throw std::runtime_error("unknown taxon '" + t.label[v] + "' in gene tree " + std::to_string(g + 1));
                int& slot = f.leaf_of[static_cast<size_t>(g) * f.n + tax];
                if (slot >= 0) throw std::runtime_error("taxon '" + t.label[v] + "' appears twice in gene tree " + std::to_string(g + 1));
                slot = off + v;
            }
            f.taxon.push_back(tax);
        }
        if (t.root != m - 1) throw std::runtime_error("internal: root is not the last node");
        f.tree_off.push_back(off + m);
    }
    return f;
}

}  // namespace wqfm
