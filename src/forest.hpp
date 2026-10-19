#pragma once

#include <algorithm>
#include <cmath>
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
    std::vector<double> supp;     // support of the branch to the parent (NaN: none given, or the root)
    std::vector<double> blen;     // length of the branch to the parent (NaN: none given, or the root)

    int nodes() const { return static_cast<int>(parent.size()); }
    int leaf(int g, int t) const { return leaf_of[static_cast<size_t>(g) * n + t]; }
};

// Two branches that are one branch of the unrooted tree (a chain through a single-child node, or
// the two branches of a degree-2 root): the larger support given, the sum of the lengths given.
inline double combine_support(double x, double y) { return std::isnan(x) ? y : std::isnan(y) ? x : std::max(x, y); }
inline double combine_length(double x, double y) { return std::isnan(x) ? y : std::isnan(y) ? x : x + y; }

// Throws std::runtime_error for unknown labels or a taxon repeated within one tree.
Forest make_forest(const std::vector<ParsedTree>& trees, const Taxa& taxa);

// Resolves every node with more than two children (the root included) arbitrarily, as the Java
// pipeline did (runResolvingPolytomy.sh: DendroPy's resolve_polytomies() without a random generator):
// the first two children in Newick order are joined under a new node, which becomes the last child,
// until two children remain. The result equals make_forest on the resolved Newick trees. New branches
// get support 0 and length 0.
Forest resolve_polytomies(const Forest& f);

}  // namespace wqfm
