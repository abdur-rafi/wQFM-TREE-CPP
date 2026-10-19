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

The profiles showed that ~90% of the time is in the per-tree kernel, mostly in subproblems with dummy
taxa. The changes, in order:

- **Internal-node and real-leaf index lists per tree.** The kernels no longer test every node for
  being a leaf. In post-order, leaves and internal nodes alternate irregularly, so that branch was
  costly.
- **Plain kernel.** Trees without dummy-taxon leaves (every tree at the top level) get a kernel
  without dummy terms.
- **Cache-line padded per-thread scratch.** Small per-thread buffers shared cache lines (false
  sharing): 8-thread gains with 4 dummies scaled 2.6× before, 3.7× after.
- **Sparse polytomy evaluator (`PolyNode`).** Polytomy nodes in trees with dummy taxa rebuilt dense
  power sums over all components for every real and dummy move, O(D·d²) each. They now keep the
  sums per component and per dummy pair, and update only what a move touches.
- **Sparse dummy taxa.** On low-resolution, high-ILS data the recursion peels off a few taxa at a
  time, so the large subproblems carry up to ~30 dummy taxa. Every dummy lies in the branch above
  every node, so each node cost O(d), or O(d²) at polytomies.
  - **Storage.** Each node now lists only the dummies with weight below it (with their weights below
    its two children, at binary nodes), replacing dense nodes × d weights.
  - **Absent dummies.** A dummy absent below a node has all its weight X in the top branch. Switching
    its side changes f there by α·X (A→B) or β·X + γ·X² (B→A):
    - α and β + γ are the top branch's one-real-taxon deltas;
    - γ sums PA over pairs of the node's other branches.

    These coefficients are summed over each tree once. At the nodes above a dummy's leaves, the
    polynomial is replaced by the exact change.
  - **Polytomy nodes** get only the listed dummies, plus the absent dummies' B-side squares on the top
    branch. Components with weight in a single branch contribute nothing else.
  - **Cost.** O(nodes + Σ list lengths) per tree instead of O(nodes · d).
- **Dummy pairs at polytomies.** Two dummies with weight in the same branch of a polytomy add a pair
  term. `PolyNode` summed these over all pairs of the p dummies listed at the node: O(p²) per node
  per pass.
  - **Pairs are dense.** On replicate 12, 98.5% of the possible pairs share a branch, because the
    top branch holds nearly every listed dummy. That amounted to ~970 million pair terms in one run.
    Storing only the pairs that occur was tried first: 18% slower at 1 thread.
  - **Sums per side instead.** Each dummy's pair terms summed over the dummies of the other side
    split into sums over its own branches of per-branch side totals (Σ T_l x_li and Σ x_li²). The
    exception is the squared term Σ_l m1_kl², which needs the side's Gram sums
    G_ij = Σ_l x_li x_lj. Their diagonal is again a per-branch total, and the others are kept only
    for the pairs of branches one dummy spans (sorted, looked up by binary search). A dummy on the
    summed side is in these totals itself, so its term with itself is subtracted.
  - **Cost.** O(entries + spanned branch pairs) per node instead of O(p²). The formula is in
    `PolyNode::pair_terms`.
- **Candidate splits scored in one pass.** The initial bipartition scored each candidate split of
  the consensus tree with its own parallel pass over the gene trees. The consensus now lists the
  distinct candidates first, each as a fresh partition or as the real taxa that move from the
  previous one (the same choice as before), and `Scorer::scores` runs every chunk of trees through
  the whole list in one parallel loop, while the chunk is in cache. The kernels take the partition
  as an argument for this, since chunks work on different candidates at the same time.

Verification:
- **Bit-identical steps.** The index lists, plain kernel and padded scratch were checked to give
  bit-identical results (`tests/bench_scorer.cpp` prints a checksum of every score and gain).
- **The two sparse changes** (polytomy evaluator, sparse dummies) change rounding only. They were
  checked by brute force (6 × 1,000 random subproblems and their children), by mutation tests (8
  deliberate bugs, all caught), and by replaying all 14 Java traces.
- **One-pass candidates are bit-identical.** Every sum is formed in the same order as before. On
  their own, they give full traces (every candidate score, gain and move) byte-identical to the
  previous build on replicate 12 (200 × 200), on 1000 × 200 and on the plant data. `brute_check` now
  also scores random batches of fresh partitions and real-taxon moves: equal bit for bit to the same
  steps through `set_partition` / `swap_real` / `score`, with the same state afterwards. It caught
  all 3 deliberate bugs planted for the test.
- **Pair terms per side change rounding only.**
  - Brute force: 6 × 1,000 random subproblems and their children.
  - Deliberate bugs: 7, all caught (4 in the first version, 3 in the final one).
  - The 10 Java trace replays on the SSD all reproduce Java's tree.
  - Against the previous build on the same three datasets: the same trees and the same move at
    every FM step, every candidate score and gain within 8e-15 relative.
  - Against the sweep: all eight conditions, replicates 1–5, 50 and 200 genes, plus replicates 1–2
    of both 500K conditions with 1000 genes. All 84 trees are byte-identical to the sweep's
    (`runs/astral2-final-check` on the SSD).
- **Against the sweep's build.** Replicates 1–5 of every condition and gene count (120 runs) give
  the same trees.
  - One differs only in Newick layout. An FM pass's end was decided at the 1e-5 boundary by
    rounding (a full swap), which relabels the two sides but leaves the tree unchanged.

| Kernel, 1000 taxa × 200 genes (ns per node) | 1 thread: score / gains | 8 threads: score / gains |
|---|---|---|
| Top level, start of session | 3.58 / 8.30 | 0.70 / 1.60 |
| Top level, now | 0.76 / 4.29 | 0.22 / 0.97 |
| 31 dummies (100 real taxa, 30 small dummies and one large): before the sparse dummies | 16.70 / 51.25 | 4.03 / 10.25 |
| 31 dummies, now | 8.96 / 25.73 | 2.01 / 5.42 |
| 1 or 4 large dummies, now | 2.5–4.7 / 7.7–12.3 | 0.6–1.1 / 1.6–2.7 |

The sparse dummies cost 10–80% on synthetic subproblems with 1–4 dummies of hundreds of taxa each, but
that does not show in whole runs:

| Whole run, same tree in every column | Start of session | Commit 1bdf66b | Now |
|---|---|---|---|
| 1000 × 200, replicate 01, 1 thread | 56.4 s | 27.8 s | 23.0 s |
| 1000 × 200, replicate 01, 8 threads | 14.2 s | 8.2 s | 7.5 s |
| 1000 × 1000, replicate 01, 8 threads (Java: 12,275 s, paper) | 167 s | 88.0 s | 88.4 s |
| 200 × 200, 500K/1e-6 replicate 12 (up to 33 dummies), 1 thread (Java: 451 s) | 71 s | 34.9 s | 16.6 s |
| 200 × 1000, 500K/1e-7 replicate 20 (the slowest of the sweep), 8 threads | – | 381 s | 165 s |
| 1KP plant, 8 threads (published tree each time) | 42.4 s | 17.0 s | 18.0 s |

Peak memory on replicate 12 halved (171 → 82 MB at 8 threads). Times vary by about ±10% between
repeated runs. Interleaved runs of both builds on high-ILS replicates 01–05 (1000 genes) took 0.8–1.0×
the earlier time. The sweep in section 3 used the build of commit 476ea53 (same results and speed as 1bdf66b); its high-ILS times are now roughly halved.

**Oversubscription.** Every FM move ends with a parallel loop over the gene trees. When more
threads than free cores are running, a descheduled thread stalls the whole move, and runs slow down
by ~10× rather than 2×. Run concurrent jobs with `-t` set so the total does not exceed the core
count.

### Thread scheduling

After the sparse dummies, 8 threads ran replicate 12 (200 × 200) only 1.8× faster than 1 thread, and
2 threads 1.25×. Profiles of multi-threaded runs showed two causes:

- **Waiting threads slept through parallel loops.** After a divide step, a thread that finishes one
  half waits for the other half, which runs as a task on another thread. A new parallel loop woke
  only the idle workers, not these waiting threads. So with 2 threads the second thread slept while
  the other half ran loop after loop: the main thread waited half of the time. `parallel_for` now
  wakes the waiting threads too (`submit` already did).
- **Moves were applied on one thread.** An FM move updated the totals along the moved taxon's path
  in every gene tree on the calling thread, with the other threads idle. At the top level of
  1000 × 200 that was about a quarter of the main thread's time. A move is now only recorded. The
  next pass over the trees applies it chunk by chunk, just before scoring the chunk, so the update
  runs in parallel with no extra barrier. Every total receives the same additions in the same order,
  so results are bit-identical.

Interleaved runs, each after a 30 s pause, same tree in every run:

| 8 threads unless noted | Wake waiters only | Plus deferred moves |
|---|---|---|
| 200 × 200, replicate 12, 4 threads | 5.68 / 5.70 s | 5.42 / 5.43 s |
| 200 × 200, replicate 12 | 5.15 / 5.08 s | 4.76 / 4.76 s |
| 1000 × 200, replicate 01 | 7.77 / 7.13 s | 5.72 / 5.69 s |

Both fixes against the build before them (that build reconstructed from the same sources), 8 threads,
45 s pause before each run, same tree as each other and as the runs in the table above:

| Whole run, 8 threads | Before | After | Peak memory before / after |
|---|---|---|---|
| 1KP plant | 15.7 s | 14.5 s | 251 / 266 MB |
| 1000 × 1000, replicate 01 | 51.0 s | 41.1 s | 772 / 920 MB |
| 200 × 1000, 500K/1e-7 replicate 20 | 166.5 s | 64.5 s | 339 / 415 MB |

Peak memory depends on how many subproblems are alive at once, which depends on the scheduling. The
same build peaked at 531 MB and 772 MB on 1000 × 1000 in two runs.

The pauses matter on the test machine, a fanless MacBook Air (M3, 4 performance and 4 efficiency
cores): back-to-back multi-threaded runs throttle, and repeated runs can differ by 2× or more. The
whole-run table above was measured back to back: its 1000 × 1000 build took 88 s there and 51 s
here.

### One-pass candidates and per-side pair terms

Against the build with both scheduling fixes. Each cell lists the runs of one build in order: two
rounds, the second in reverse build order, 30 s pause before every run, same tree in every run:

| Whole run | Before | One-pass candidates | Plus per-side pair terms |
|---|---|---|---|
| 200 × 200, replicate 12, 1 thread | 16.55 / 16.46 / 16.61 / 16.50 s | 16.21 / 16.16 s | 14.38 / 14.38 s |
| 200 × 200, replicate 12, 8 threads | 4.62 / 4.63 s | 4.18 / 4.15 s | 3.79 / 3.85 s |
| 1000 × 200, replicate 01, 1 thread | 22.48 / 22.27 s | – | 22.01 / 21.94 s |
| 1000 × 200, replicate 01, 8 threads | 5.41 / 5.41 s | 5.06 / 4.95 s | 4.97 / 4.99 s |
| 1KP plant, 8 threads | 14.38 / 14.56 / 17.55 / 15.56 s | 14.79 / 14.88 s | 14.61 / 14.59 s |
| 200 × 1000, 500K/1e-7 replicate 20, 8 threads (one round) | 52.75 s | 51.19 s | 48.91 s |
| 1000 × 1000, replicate 01, 8 threads (one round) | 40.88 s | – | 35.90 s |

The last two rows used the first version of the pair terms. That version cleared per-branch arrays at
every polytomy node, and on plant (polytomies of up to 576 branches, few dummies) it took
15.25 / 15.28 s, slower than one-pass candidates alone. The final version works over the entries
only. Its results are bit-identical to the first version, and it is faster on both plant and
replicate 12. One-pass candidates matter most with threads, where they replace a parallel loop per
candidate. The pair terms matter where high-ILS polytomies carry many dummies.

### Memory

Peak memory grew with the thread count by up to 60% (1000 × 1000: 557 MB with 1 thread, 781–908 MB
with 8) and varied between runs. Two causes:

- **Reused scorers kept their largest buffers.** Each subproblem has a scorer holding its compressed
  forest. Scorers were pooled and reused, and a reused scorer kept the capacity of the largest forest
  it had ever held (often the top level's). More threads put more scorers in circulation, and which
  one went where depended on scheduling. Scorers are now created per subproblem and freed when it is
  done; that costs no measurable time.
- **The parsed gene trees stayed in memory.** The Newick parser's per-node vectors and label strings
  (~60 bytes a node, ~120 MB at 1000 × 1000) were kept for the whole run, although only the compact
  forest is used after setup. They are now freed once the forest is built.

Peak memory footprint (macOS), before / after, two runs each, same tree in every run:

| Run | Before | After |
|---|---|---|
| 1000 × 200, 1 thread | 116 / 115 MB | 84 / 89 MB |
| 1000 × 200, 8 threads | 133 / 120 MB | 95 / 104 MB |
| 1KP plant, 8 threads | 300 / 297 MB | 196 / 165 MB |
| 1000 × 1000, 1 thread | 510 / 532 MB | 285 / 295 MB |
| 1000 × 1000, 8 threads | 572 / 620 MB | 361 / 324 MB |

With 8 threads the peak is now 15–20% above the single-threaded one. Measure memory with the peak
footprint, not the maximum resident set: under memory pressure or when the laptop sleeps, macOS
compresses or swaps pages and the resident set understates the peak (the same 1-thread run showed
303 MB and 406 MB resident).
