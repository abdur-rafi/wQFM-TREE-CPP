#pragma once

#include <string>

#include "forest.hpp"

namespace wqfm {

// Greedy (majority-rule extended) consensus of the gene trees, the role PAUP's
// "contree all / majrule le50" plays in the Java pipeline (run.sh). Splits are added in order of
// decreasing frequency (ties: first seen, in gene-tree order), each kept if compatible with all
// accepted splits, until the tree is fully resolved.
//
// PAUP may break frequency ties differently, which can change the initial bipartition of a
// subproblem in rare cases; pass PAUP's consensus with -c to reproduce a Java run exactly.
// Requires complete gene trees (every taxon in every tree); throws otherwise.
std::string greedy_consensus_newick(const Forest& f, const Taxa& taxa);

}  // namespace wqfm
