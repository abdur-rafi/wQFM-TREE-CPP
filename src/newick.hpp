#pragma once

#include <string>
#include <vector>

namespace wqfm {

// A parsed tree with nodes in creation order: a leaf is created when its label is read and an
// internal node at its closing ')', so children always precede their parent and the root is last.
// Children are kept in REVERSE input order, as in the Java implementation (its parser pops them off
// a stack). This matters only where the Java code iterates children: the consensus-candidate order.
struct ParsedTree {
    std::vector<int> parent;                  // -1 for the root
    std::vector<std::vector<int>> children;
    std::vector<std::string> label;           // leaf labels; empty for internal nodes
    int root = -1;
};

// Accepts branch lengths, internal labels / support values, [comments] and 'quoted labels'.
// Throws std::runtime_error on malformed input.
ParsedTree parse_newick(const std::string& text);

}  // namespace wqfm
