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
  - **The C++ code does the same at larger scales.** In the default run on 1000 taxa × 1000 genes
    (replicate 01), 13 FM passes kept the prefix that moves every unit, with cumulative gains of
    about 5e-4. That is the same bipartition, so the true gain is 0, and each such pass is wasted.
  - **`--relative-eps` (opt-in) avoids it.** FM then keeps a longer prefix only if it beats the
    best by 1e-5 + 1e-12 × (sum of |gain| of the pass so far).
    - Test: 8 replicates across the conditions (50 to 1000 genes, 200 to 1000 taxa).
    - Result: the same trees as the default. On the two 1000 × 1000 replicates, 30% and 6% fewer
      FM moves (19,981 → 13,952 and 25,065 → 23,623).

### Consensus

The published trees used PAUP's greedy consensus; the C++ runs use the built-in one. Where they
differ, so can the species tree.
- **Replicate 11 (500K, 1e-6, 1000 genes).** Java and C++ with the built-in consensus give the same
  tree, 7 splits from the published one.
- **Replicate 03 (50 genes).** Java with the built-in consensus is 16 splits from the published tree.

With 50 genes and high ILS, split frequencies tie often, and PAUP breaks those ties in its own order.

## 2. Accuracy against the true trees

All replicates with a published tree (1,101 runs), built-in consensus, 8 threads. Each cell: mean
missing-branch rate against the true tree for the C++ run / for the published wQFM-TREE tree, and
in brackets the number of replicates where the two trees are identical.

| Condition (taxa, height, rate) | Reps | 50 genes | 200 genes | 1000 genes |
|---|---|---|---|---|
| 200 taxa, 500K, 1e-6 | 47 | 0.1879 / 0.1892 (1) | 0.0952 / 0.0971 (6) | 0.0537 / 0.0542 (11) |
| 200 taxa, 500K, 1e-7 | 50 | 0.2052 / 0.2041 (2) | 0.1088 / 0.1071 (1) | 0.0653 / 0.0632 (5) |
| 200 taxa, 2M, 1e-6 | 50 | 0.0801 / 0.0793 (33) | 0.0471 / 0.0467 (43) | 0.0317 / 0.0316 (47) |
| 200 taxa, 2M, 1e-7 | 50 | 0.0968 / 0.0973 (15) | 0.0575 / 0.0579 (18) | 0.0418 / 0.0410 (26) |
| 200 taxa, 10M, 1e-6 | 50 | 0.0499 / 0.0500 (41) | 0.0273 / 0.0273 (50) | 0.0174 / 0.0175 (49) |
| 200 taxa, 10M, 1e-7 | 50 | 0.0675 / 0.0677 (39) | 0.0490 / 0.0492 (42) | 0.0425 / 0.0428 (44) |
| 500 taxa, 2M, 1e-6 | 50 | 0.0752 / 0.0754 (15) | 0.0430 / 0.0431 (28) | 0.0284 / 0.0280 (26) |
| 1000 taxa, 2M, 1e-6 | 20 | 0.0808 / 0.0808 (0) | 0.0480 / 0.0477 (3) | 0.0327 / 0.0325 (13) |

- **Same accuracy as published.** Over all runs, the C++ tree has a missing-branch rate +0.00009
  higher on average: better in 191 runs, worse in 217, equal in 693, and identical to the published
  tree in 558.
- **Where the trees differ, the consensus is the cause** (section 1). Identical trees are rarest
  where frequency ties are commonest: high ILS (500K) and 1000 taxa with 50 genes.

## 3. Speed

Same sweep (8 threads, Apple M3). The Java times and memory are Supplementary Tables 4 and 5:
- **Different setup.** An Intel i7-10700K, averaged over 10 replicates (200 and 500 taxa) or 5
  (1000 taxa).
- **200-taxon condition assumed.** Table 4 does not name it; it is assumed to be the default 2M, 1e-6.

| Condition | Genes | C++ 8 threads: mean (max) | Java, paper | Ratio of means | C++ peak memory (max) | Java memory, paper |
|---|---|---|---|---|---|---|
| 200 taxa, 500K, 1e-6 | 50 | 1.84 s (14.3) | – | – | 67 MB | – |
| 200 taxa, 500K, 1e-6 | 200 | 4.98 s (23.4) | – | – | 194 MB | – |
| 200 taxa, 500K, 1e-6 | 1000 | 26.65 s (129.8) | – | – | 432 MB | – |
| 200 taxa, 500K, 1e-7 | 50 | 2.97 s (18.9) | – | – | 79 MB | – |
| 200 taxa, 500K, 1e-7 | 200 | 9.30 s (61.9) | – | – | 225 MB | – |
| 200 taxa, 500K, 1e-7 | 1000 | 50.23 s (330.3) | – | – | 481 MB | – |
| 200 taxa, 2M, 1e-6 | 50 | 0.18 s (0.3) | 20.1 s | 113× | 18 MB | – |
| 200 taxa, 2M, 1e-6 | 200 | 0.52 s (0.9) | 93.9 s | 181× | 61 MB | – |
| 200 taxa, 2M, 1e-6 | 1000 | 2.38 s (4.3) | 439.5 s | 184× | 278 MB | 2.0–2.6 GB |
| 200 taxa, 2M, 1e-7 | 50 | 0.33 s (2.5) | – | – | 26 MB | – |
| 200 taxa, 2M, 1e-7 | 200 | 1.00 s (7.3) | – | – | 145 MB | – |
| 200 taxa, 2M, 1e-7 | 1000 | 4.94 s (40.1) | – | – | 393 MB | – |
| 200 taxa, 10M, 1e-6 | 50 | 0.18 s (0.5) | – | – | 16 MB | – |
| 200 taxa, 10M, 1e-6 | 200 | 0.53 s (1.1) | – | – | 48 MB | – |
| 200 taxa, 10M, 1e-6 | 1000 | 2.52 s (5.2) | – | – | 187 MB | – |
| 200 taxa, 10M, 1e-7 | 50 | 0.22 s (1.4) | – | – | 22 MB | – |
| 200 taxa, 10M, 1e-7 | 200 | 0.62 s (3.1) | – | – | 74 MB | – |
| 200 taxa, 10M, 1e-7 | 1000 | 2.97 s (17.7) | – | – | 303 MB | – |
| 500 taxa, 2M, 1e-6 | 50 | 1.04 s (2.6) | 156.6 s | 150× | 52 MB | – |
| 500 taxa, 2M, 1e-6 | 200 | 3.27 s (8.2) | 793.5 s | 243× | 171 MB | – |
| 500 taxa, 2M, 1e-6 | 1000 | 19.52 s (45.4) | 3,413.7 s | 175× | 622 MB | 4.9 GB |
| 1000 taxa, 2M, 1e-6 | 50 | 4.30 s (15.1) | 943.2 s | 219× | 166 MB | – |
| 1000 taxa, 2M, 1e-6 | 200 | 16.01 s (49.2) | 3,155.8 s | 197× | 379 MB | – |
| 1000 taxa, 2M, 1e-6 | 1000 | 86.30 s (271.6) | 12,275.2 s | 142× | 1,124 MB | 6.7 GB |

- **The slowest runs are the high-ILS conditions (500K).** Their gene trees are poorly resolved, and
  the recursion builds subproblems with many dummy taxa (see "Where the time still goes").
- **Peak memory grows with the number of dummy taxa.** The dense per-node dummy weights take
  nodes × d doubles.

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
| Top level: score / gains (ns per node) | 3.58 → 0.77 / 8.30 → 4.33 | 0.70 → 0.23 / 1.60 → 0.94 |
| 1 dummy: score / gains | 4.17 → 2.09 / 10.68 → 6.97 | 0.89 → 0.50 / 2.42 → 1.44 |
| 4 dummies: score / gains | 6.87 → 2.49 / 19.23 → 9.23 | 1.45 → 0.65 / 4.30 → 1.91 |

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
