#pragma once

#include <cstdint>

#include "forest.hpp"

namespace wqfm {

// Quartet support of a species tree (ASTRAL's quartet score), counted exactly: over all gene trees,
// the quartets a gene tree resolves (on its own taxa), and how many of them the species tree resolves
// the same way. Trees are unrooted; polytomies and missing taxa are allowed.
struct QuartetCounts {
    int64_t agree = 0;      // gene-tree quartets that the species tree shows the same way
    int64_t resolved = 0;   // quartets resolved by the gene trees
};

// `species` holds one tree, containing every taxon of `genes`. O(nodes(species) * nodes(gene tree))
// per gene tree; gene trees are split over `threads`.
QuartetCounts quartet_support(const Forest& genes, const Forest& species, int threads = 1);

}  // namespace wqfm
