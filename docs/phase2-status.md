# Phase 2 status

Single-threaded optimizations of the Phase 1 search. Each one was checked by brute force (with
planted-bug tests), by replaying Java traces, and by independent runs against Java.

## Changes

1. **Compressed gene trees per subproblem** (`SubForest`, `Scorer::compress_from`). Every maximal part
   of a gene tree whose taxa all lie under one dummy taxon becomes one weighted leaf. Purity is judged
   on the unrooted tree, from below and above each node, and same-dummy siblings are merged. Gene
   trees with fewer than 4 units are dropped. This is exact (see `src/scorer.hpp`).
2. **Child forests derived from the parent's** (`Scorer::setup_child`). The child's units are unions
   of the parent's, so compressing the parent's compressed forest is exact and costs O(parent size)
   instead of O(all gene-tree nodes). Setup time on the 1KP data fell from 15.8 s to 0.8 s.
3. **One setup per subproblem.** Each subproblem's scorer is built once and used first for its
   initial bipartition, then for its FM passes.
4. **Incremental scoring of initial-split candidates.** Consecutive consensus candidates usually
   differ in a few real taxa, so the scorer moves those taxa instead of rebuilding its state.
5. **Fast path for polytomy nodes in gene trees without dummy taxa** (the whole top level). It uses
   fixed-size power sums, the formula that was bit-exact against Java in Phase 0.
6. **Single-child nodes are contracted when reading gene trees.** Topology is unchanged; the Java
   version exits on such nodes.

`--no-compress` switches 1–2 off, for testing.

## Validation

- **`brute_check`** compares the compressed path, the plain path and derived child forests
  (~590 per run) against quartet enumeration: scores agree to ~1e-13 and gains to ~1e-12. Planted
  bugs in cluster weights and in the derived-child weights are caught.
- **Replay of Java traces** on all four simulated datasets: same trees, every candidate score and
  gain within 1e-11, and 0 real disagreements.
- **Independent runs** with Java's taxon order give Java's tree on all four datasets. On the 1KP
  data the result is the published wQFM-TREE tree (`plant/plant.ewqfm.A.tre`, RF = 0).

## Timings

Single-threaded on an Apple M3, run one at a time. A long Java analysis was running on another
core. Java times for the simulated sets come from the traced validation runs.

| Dataset | Taxa × genes | Java | Phase 1 | Phase 2 | Phase 2 vs Java |
|---|---|---|---|---|---|
| Simulated, discordance 0.44 | 200 × 200 | 44.6 s | 6.78 s | 1.07 s | 42× |
| Simulated, discordance 0.77 | 200 × 200 | 82.6 s | 7.30 s | 1.71 s | 48× |
| 15% missing taxa + polytomies | 120 × 150 | 25.9 s | 2.59 s | 0.51 s | 51× |
| Simulated, discordance 0.34 | 300 × 150 | 84.6 s | 10.02 s | 1.61 s | 53× |
| 1KP plant | 1,178 × 410 | ≈5.5 h (paper) | 1,203 s | 174 s | ≈110× |

On the 1KP data the 174 s break down as initial bipartitions 34.9 s, FM passes 138.1 s and setup
0.8 s. Peak memory is 124 MB; the Java version needs ~5.7 GB (Supp. Table 5).

## Measured and deferred

- **Sparse dummy-taxon gains.** Node visits see on average only 0.75–1.0 present dummy taxa (and
  0.1–0.3 below the node), because the top levels dominate. Under replay, dummy gains are 5–13% of
  FM time, so the sparse formulation (whose top-only part is linear/quadratic in the dummy's
  weight, see `docs/plan.md`) was not worth its complexity yet.
- **Kernel speed.** About 14 ns per node visit single-threaded, close to the Phase 0 prototype.
  Further single-thread gains need SIMD batching or a different data layout; parallelism (Phase 3)
  is the larger lever.
- **Integer top level.** Not needed: at the top level every value is an integer below 2^53, so
  doubles are already exact.
- **Duplicate compressed trees** (stored once with a count) and **batched candidate scoring**:
  candidates for later.
