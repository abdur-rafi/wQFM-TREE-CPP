#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "newick.hpp"

namespace wqfm {

struct Taxa {
    std::vector<std::string> label;
    std::unordered_map<std::string, int> index;

    int add(const std::string& l);           // id of l, adding it if new
    int find(const std::string& l) const;    // -1 if unknown
    int size() const { return static_cast<int>(label.size()); }
};

// All gene trees in one flat array. Within a tree, children precede their parent and the root is the
// last node; a degree-2 root is fine (it contributes nothing and its children see each other as "top").
struct Forest {
    int n = 0;                    // taxa
    int k = 0;                    // trees
    std::vector<int> tree_off;    // k+1 node offsets
    std::vector<int> parent;      // -1 for a tree root
    std::vector<int> child_off;   // nodes+1, CSR offsets into children
    std::vector<int> children;
    std::vector<int> taxon;       // taxon id of a leaf, -1 for internal nodes
    std::vector<int> leaf_of;     // k*n: node of taxon t in tree g, or -1 if absent

    int nodes() const { return static_cast<int>(parent.size()); }
    int leaf(int g, int t) const { return leaf_of[static_cast<size_t>(g) * n + t]; }
};

// Throws std::runtime_error for unknown labels or a taxon repeated within one tree.
Forest make_forest(const std::vector<ParsedTree>& trees, const Taxa& taxa);

}  // namespace wqfm
