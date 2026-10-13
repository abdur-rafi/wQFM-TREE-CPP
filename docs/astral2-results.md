# ASTRAL-II datasets: validation, accuracy and speed

The paper's simulated benchmark (Mirarab & Warnow 2015): the eight model conditions analysed in the
paper, 50, 200 and 1000 estimated gene trees, against the true species trees and the published
wQFM-TREE trees in `../results/Astral_2_datasets`.

Data: Dryad `estimated-gene-trees.tar.bz`, unpacked outside the repository
(`<data>/<condition>/<replicate>/estimatedgenetre`, 1000 trees per file). The 50- and 200-gene
analyses use the first 50 and 200 trees of each file; the published trees confirm this (most are
reproduced exactly).

Scripts:
- `tools/astral2/run.py`: runs the C++ version over conditions, replicates and gene counts, and
  scores each tree against the true tree and the published one.
- `tools/astral2/validate_batch.sh`: runs `tools/compare/validate.sh` (traced Java, independent C++
  run, C++ replay of Java's decisions) on chosen replicates, with the same gene trees and consensus
  for both.

Machine: Apple M3 (4 performance + 4 efficiency cores), 16 GB.

## 1. Same results as Java

| Condition, replicate 01 | Genes | Java (traced) | C++ (traced) | Replay | Independent C++ run, same consensus and taxon order |
|---|---|---|---|---|---|
| 200 taxa, 500K, 1e-6 | 50 | 12.1 s | 0.6 s | agrees | identical tree |
| 200 taxa, 10M, 1e-6 | 50 | 4.8 s | 0.2 s | agrees | 1 split (Java's rounding: see below) |
| 200 taxa, 2M, 1e-7 | 50 | 7.3 s | 0.3 s | agrees | 1 split (FM near-tie) |
| 200 taxa, 500K, 1e-7 | 50 | 15.5 s | 0.9 s | agrees | identical tree |
| 500 taxa, 2M, 1e-6 | 50 | 82.7 s | 2.0 s | agrees | 1 split (FM near-tie) |
| 200 taxa, 500K, 1e-6 | 200 | 87.8 s | 2.3 s | agrees | identical tree |
| 200 taxa, 500K, 1e-7 | 200 | 69.4 s | 2.1 s | agrees | identical tree |
| 200 taxa, 500K, 1e-6 | 1000 | 388.8 s | 6.9 s | agrees | identical tree |
| 1000 taxa, 2M, 1e-6 | 50 | 412.9 s | 8.6 s | agrees | identical tree |

- **Replay** means the C++ run follows Java's decisions and logs its own numbers. In every case it
  reproduces Java's tree, all 74,109 candidate scores and 45,666 FM gains agree within 8e-12
  (relative to the largest value in the candidate set or FM pass), and the C++ code never prefers
  another move except at exact ties.
- **Replicate 03** of the highest-ILS condition (50 genes; Java 77 s) also agrees under replay.

### Ties

- **Initial bipartition.** In 15–45% of subproblems several consensus candidates have the same
  score. Java's tied scores are bit-identical (same quartet weights summed), so its strict `>`
  keeps the first. The C++ scores carried rounding noise, which picked an arbitrary one.
  `ConsensusPartitioner::make` now treats scores within `1e-9 · max|score| + 1e-16 · Σ_g W_g⁴` of the
  best as tied and keeps the first. On all traces, noisy ties differ by at most 4e-14 of the largest
  |score|, and distinct best scores by at least 1.4e-4.
  - Effect: independent runs now match Java's initial split in every subproblem (up to 18 per run
    differed before), and three of the eight runs above became identical to Java's tree.
- **Remaining differences are Java's own rounding.** Examples:
  - FM moves with equal gains. Java picked the earlier move 11 times and the later 12 times,
    depending on its rounding.
  - A candidate set where every score is 0 but one came out as +1.4e-14 in Java.
  - Java's EPS = 1e-5 for keeping an FM prefix is absolute. With cumulative gains of order 1e7–1e8,
    Java's rounding exceeds it: a full pass that moves every unit (the same bipartition, true gain
    0) came out at +2.9e-4 and counted as an improvement, starting further passes. The C++ sum was
    3.3e-8.

### Consensus

The published trees used PAUP's greedy consensus; the C++ runs use the built-in one. Where they
differ, so can the species tree.
- **Replicate 11 (500K, 1e-6, 1000 genes).** Java and C++ with the built-in consensus give the same
  tree, 7 splits from the published one.
- **Replicate 03 (50 genes).** Java with the built-in consensus is 16 splits from the published tree.

With 50 genes and high ILS, split frequencies tie often, and PAUP breaks those ties in its own order.

## 2. Accuracy against the true trees

*Filled in from `runs/astral2-cpp-t8/summary.txt` (8 threads, all published replicates).*

## 3. Speed

*Per-condition times from the same sweep; Java times from Supplementary Table 4 (Intel i7-10700K).*

### Kernel work during this analysis

All of it was checked to give bit-identical results (`tests/bench_scorer.cpp` prints a checksum of every
score and gain), except the polytomy evaluator, which was checked by brute force and Java replays.
The profiles showed that ~90% of the time is in the per-tree kernel, mostly in subproblems with dummy
taxa.

- **Internal-node and real-leaf index lists per tree.** The kernels no longer test every node for
  being a leaf. In post-order, leaves and internal nodes alternate irregularly, so that branch was
  costly.
- **Plain kernel.** Trees without dummy-taxon leaves (every tree at the top level) get a kernel
  without dummy terms.
- **Kernel specialized for 1–4 present dummies.** Per-dummy values and gains stay in registers.
- **Cache-line padded per-thread scratch.** Small per-thread buffers shared cache lines (false
  sharing): 8-thread gains with 4 dummies scaled 2.6× before, 3.7× after.
- **Sparse polytomy evaluator (`PolyNode`).** Polytomy nodes in trees with dummy taxa rebuilt dense
  power sums over all components for every real and dummy move, O(D·d²) each. They now keep the
  sums per component and per dummy pair and update only what a move touches.

| Kernel, 1000 taxa × 200 genes | 1 thread, before → after | 8 threads, before → after |
|---|---|---|
| Top level: score / gains (ns per node) | 3.58 → 0.76 / 8.30 → 4.29 | 0.70 → 0.22 / 1.60 → 0.96 |
| 1 dummy: score / gains | 4.17 → 2.10 / 10.68 → 7.16 | 0.89 → … / 2.42 → … |
| 4 dummies: score / gains | 6.87 → 2.47 / 19.23 → 9.46 | 1.45 → … / 4.30 → … |

| Whole run | Before | After |
|---|---|---|
| 1000 × 200, replicate 01, 1 thread | 56.4 s | 28.8 s |
| 1000 × 200, replicate 01, 8 threads | 14.2 s | 9.1 s |
| 1000 × 1000, replicate 01, 8 threads | 167 s | 79.6 s (Java: 12,275 s, paper) |
| 200 × 200, 500K/1e-6 replicate 12 (8,785 polytomies, up to 33 dummy taxa), 1 thread | 71 s | 34.8 s (Java: 451 s) |

All runs gave the same trees before and after.

**Where the time still goes.** On low-resolution, high-ILS data the recursion peels off a few taxa
at a time, so the large subproblems carry up to ~30 dummy taxa. The cost per node is then O(d) for
binary nodes and O(d²) for polytomies. This is because the branch above each node (the rest of the
tree) contains nearly every dummy. A fully sparse formulation would make both O(1 + dummies below
the node):
- Gram sums over branch pairs for the dummy-pair terms.
- A per-tree polynomial in each dummy's total weight for the nodes it is absent from below.

That is not implemented.

**Oversubscription.** Every FM move ends with a parallel loop over the gene trees. When more
threads than free cores are running, a descheduled thread stalls the whole move, and runs slow down
by ~10× rather than 2×. Run concurrent jobs with `-t` set so the total does not exceed the core
count.
