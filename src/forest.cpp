#include "forest.hpp"

#include <deque>
#include <stdexcept>
#include <utility>

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
            // the branch to the parent; a chain through contracted single-child nodes is one branch
            double s = t.support[v], l = t.length[v];
            int u = t.parent[v];
            while (u >= 0 && t.children[u].size() == 1) {
                s = combine_support(s, t.support[u]);
                l = combine_length(l, t.length[u]);
                u = t.parent[u];
            }
            f.supp.push_back(u < 0 ? std::nan("") : s);
            f.blen.push_back(u < 0 ? std::nan("") : l);
        }
        f.tree_off.push_back(off + kept);
    }
    return f;
}

Forest resolve_polytomies(const Forest& f) {
    Forest r;
    r.n = f.n;
    r.k = f.k;
    r.tree_off.push_back(0);
    r.child_off.push_back(0);
    r.leaf_of.assign(f.leaf_of.size(), -1);
    std::vector<std::vector<int>> kids;     // children in Newick order; ids >= m are new nodes
    std::vector<int> id;
    std::vector<std::pair<int, size_t>> st;
    std::deque<int> q;
    for (int g = 0; g < f.k; ++g) {
        const int off = f.tree_off[g], m = f.tree_off[g + 1] - off;
        kids.assign(m, {});
        for (int v = 0; v < m; ++v)          // forest children are in reverse Newick order
            for (int i = f.child_off[off + v + 1] - 1; i >= f.child_off[off + v]; --i) kids[v].push_back(f.children[i] - off);
        for (int v = 0; v < m; ++v) {
            if (kids[v].size() <= 2) continue;
            q.assign(kids[v].begin(), kids[v].end());
            while (q.size() > 2) {
                const int a = q.front(); q.pop_front();
                const int b = q.front(); q.pop_front();
                q.push_back(static_cast<int>(kids.size()));
                kids.push_back({a, b});
            }
            kids[v].assign(q.begin(), q.end());
        }
        // renumber in postorder, as the parser numbers the resolved Newick tree
        id.assign(kids.size(), -1);
        st.assign(1, {m - 1, 0});
        while (!st.empty()) {
            const int v = st.back().first;
            if (st.back().second < kids[v].size()) {
                const int c = kids[v][st.back().second++];
                st.emplace_back(c, 0);
                continue;
            }
            st.pop_back();
            const int nv = r.nodes();
            id[v] = nv;
            for (auto it = kids[v].rbegin(); it != kids[v].rend(); ++it) {
                r.children.push_back(id[*it]);
                r.parent[id[*it]] = nv;
            }
            r.child_off.push_back(static_cast<int>(r.children.size()));
            r.parent.push_back(-1);
            const int tax = v < m ? f.taxon[off + v] : -1;
            r.taxon.push_back(tax);
            r.supp.push_back(v < m ? f.supp[off + v] : 0.0);
            r.blen.push_back(v < m ? f.blen[off + v] : 0.0);
            if (tax >= 0) r.leaf_of[static_cast<size_t>(g) * r.n + tax] = nv;
        }
        r.tree_off.push_back(r.nodes());
    }
    return r;
}

}  // namespace wqfm
