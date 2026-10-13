# Phase 1 status

The complete wQFM-TREE search is implemented in C++, single-threaded: `build/wqfm-tree`.

| File | Contents |
|---|---|
| `src/newick.*` | Newick parser (branch lengths, support values, comments, quoted labels); keeps Java's reversed child order |
| `src/forest.*` | taxa, and all gene trees in one flat array |
| `src/units.hpp` | dummy taxa (nested structure, Eq. 1 weights) and subproblems |
| `src/scorer.*` | score and gains of real and dummy taxa from the gene trees; closed-form polytomy nodes (`PolyEval`) |
| `src/consensus.*` | initial bipartition from the consensus tree, replicating Java's candidate order and root re-balancing |
| `src/wqfm.*` | FM passes (Java's tie-breaking, ≤ 6 passes), divide, recursion, joining subtrees |
| `src/trace.*` | step-by-step trace (shared format with the Java build) and replay of a Java trace |

## How it is validated

1. **Brute force** (`build/brute_check`). The score, every real and dummy gain, and the state after
   incremental swaps are compared with explicit quartet enumeration. Subproblems are random, with
   polytomies, missing taxa, rooted and unrooted gene trees, nested dummy taxa, and some with no real
   taxa at all. Agreement is ~1e-13. Planted single-term bugs are all caught.
2. **Independent runs.** The traced Java build and the C++ program run on the same data with the same
   taxon order, and their traces are compared subproblem by subproblem.
3. **Replay.** The C++ program follows the decisions in a Java trace while computing and tracing its own
   values. Every candidate score and every gain at every level is therefore compared, even after a
   near-tie that rounding broke differently. Where the C++ would have chosen another move, the gap
   between the two gains separates near-ties from real disagreements.

`tools/compare/validate.sh` runs layers 2 and 3 on any dataset; `make check` includes a small instance.

## Results

Simulated data come from `tools/sim` (multispecies coalescent; discordance = mean RF between gene
and species trees). Times are for traced runs on an Apple M3, under the background load described in
[phase0-results.md](phase0-results.md).

| Dataset | Taxa × genes | Java | C++ | Independent run vs Java | Replay |
|---|---|---|---|---|---|
| 11-taxon example (repo `input/`) | 11 × 5 | – | – | 1 split differs: Java broke an exact tie (both gains −1) by rounding noise | same tree; 91 scores, 40 gains within 4.4e-15 |
| Simulated, discordance 0.44 | 200 × 200 | 44.6 s | 6.4 s | **same tree**; 1 near-tie broken differently, no effect | same tree; within 5.4e-13 |
| Simulated, discordance 0.77 | 200 × 200 | 82.6 s | 10.4 s | **same tree**, every move identical | same tree; within 8.8e-12 |
| 15% missing taxa + polytomies | 120 × 150 | 25.9 s | 3.7 s | **same tree** (a near-tie in one initial split, same final splits) | same tree; within 3.7e-12 |
| Simulated, discordance 0.34 | 300 × 150 | 84.6 s | 15.1 s | **same tree**; 3 near-ties, no effect | same tree; within 6.8e-12 |
| 1KP plant | 1,178 × 410 | in progress | in progress | | |

In every replay there were 0 real disagreements: every move the C++ would have chosen differently
was a near-tie (relative gain gap ≤ 1e-9).

## Behaviour and known differences

- **Ties.** Java breaks exact ties by rounding noise, which cannot be reproduced. The C++ computes them
  (near-)exactly and applies Java's rule: the first real taxon wins, and a dummy taxon only if strictly
  better. Trees can therefore differ when the data contain exact ties, as in the 11-taxon example.
- **Taxon order.** Taxon ids decide ties. Java numbers taxa in `HashSet` order; the C++ default is the
  order of first appearance. `--taxa-order <java trace>` reproduces Java's order.
- **Consensus.** A consensus tree is still required (`-c`), as with the Java version. The native
  greedy consensus is Phase 3.
- **Duplicate nodes.** Java scores identical gene-tree nodes once. The C++ does not deduplicate yet;
  results are the same and only speed differs (Phase 2).

## Bug found by validation

In subproblems with no real taxa, all dummy gains were skipped: "compute gains" was signalled by a
non-null `real_gain` pointer, and an empty vector's `data()` may be null. The random subproblems in
the brute-force test always kept some real taxa, so it went unnoticed until the Java comparison. The
flag is now explicit, and the test includes all-dummy subproblems (it fails on the old code).
