# Phase 0 results

**Machine:** Apple M3 (4 performance + 4 efficiency cores), 16 GB RAM, macOS; Temurin JDK 21;
Apple clang 16 (`-O3 -mcpu=native -ffp-contract=off`).

**Load during timings:** one long Java analysis was running in the background throughout (one core
busy), and the machine was under memory pressure from a virtual machine. Timings are therefore
indicative; Java and C++ were measured under the same conditions.

## 1. Correctness

| Check | Result |
|---|---|
| Java score and gains vs brute-force quartet enumeration (1,500 random subproblems: polytomies, missing taxa, nested dummies, duplicate nodes, after swaps) | agree to ~1e-13 |
| Same check with non-quartet type B | fails on every polytomy-bearing case, so **A is the paper's formula** |
| Collapsing single-dummy subtrees, scored by brute force (1,500 subproblems) | **exact** (~1e-14); `Phase0Stats` node counts match an explicit construction |
| C++ prototype vs Java at the top level | **bit-identical** scores, gains, every FM move and final partition, on all datasets below and for 1–8 threads |

## 2. Simulated data: compression and time split

200 taxa × 200 genes, Yule species tree, gene trees simulated under the multispecies coalescent
(`tools/sim`). These are true gene trees, so likely more favourable to compression than estimated ones.

| Mean gene-tree discordance (RF to species tree) | 0.06 | 0.18 | 0.28 (+ NNI error) | 0.44 | 0.77 |
|---|---|---|---|---|---|
| FM node visits, Java → compressed | 4.7× fewer | 4.8× | 4.8× | 4.9× | 4.8× |
| Initial-split node visits | 4.6× fewer | 4.6× | 4.6× | 4.6× | 4.4× |
| Share of Java runtime in the initial split | 57% | 54% | 51% | 54% | 53% |
| Candidates considered / rescored by Java / distinct | 3985/3757/3287 | 4021/3770/3299 | 4032/3790/3305 | 4080/3823/3329 | 4333/4057/3528 |
| FM passes per subproblem (mean / max) | 1.03 / 2 | 1.09 / 2 | 1.11 / 2 | 1.19 / 2 | 1.54 / 3 |
| Max dummy taxa in a subproblem (Java crash bound 62) | 5 | 6 | 6 | 6 | 6 |
| Total Java runtime | 29 s | 29 s | 39 s | 58 s | 72 s |

- **Compressed trees keep roughly s/n of the nodes** at each recursion level (s = subproblem size,
  n = taxa). At discordance 0.77 the deepest levels keep ~5.5% instead of 1–2%.
- **After compression, the top two or three levels dominate the work.** In Java, each level costs
  about the same.

## 3. 1KP plant data (410 genes × 1,178 taxa)

- **No duplicate nodes:** 233,134 unique tripartitions out of 233,134 nodes. With incomplete gene
  trees, each node's key includes that tree's present-taxon set.
- **Polytomies:** 24,417 polytomy nodes; the largest has 576 branches.
- **Top level of the Java run:** initial split 314 s (2,268 candidates, all distinct); FM 1,135 s
  over 2 passes; the two children's initial splits 365 s.
- **Compression on real data** (instrumented Java run, stopped after 940 of ~1,175 subproblems and
  about 10 hours, once the C++ version made it unnecessary): compressed trees keep roughly s/n of
  the nodes at each level (level 3: 27% at s = 296), as on simulated data. FM needs 3.9× fewer node
  visits and the initial split 5.2× fewer. Candidate dedup would save 1.1×. At most 8 dummy taxa per
  subproblem (Java's crash bound is 82).

## 4. Gain-sweep prototype

Time to compute all gains once:

| Dataset | Nodes | Java | C++, 1 thread | C++, multithreaded |
|---|---|---|---|---|
| Simulated 200 × 200 | 80K | 10.2 ms | 0.52 ms (19×) | 0.33 ms (2 threads) |
| Incomplete + polytomies, 150 × 100 | — | 5.65 ms | — | 0.081 ms (4 threads) |
| 1KP plant 1,178 × 410 | 600K | 644 ms | 8.2 ms (78×) | 1.8–2.8 ms (8 threads) |
| Simulated 1,000 × 1,000 | 2.0M | not run (memory) | 22.7 ms | 7.9 ms (4), 7.1 ms (8) |

One full FM pass at the top level:

| Dataset | Java | C++ |
|---|---|---|
| Simulated 200 × 200 | 4.10 s | 0.096 s (4 threads) |
| Incomplete + polytomies | 1.00 s | 0.015 s (4 threads) |
| 1KP plant | ≈567 s (1,135 s for 2 passes, measured in the full run) | 15.9 s (1 thread), 4.3 s (8 threads) |
| Simulated 1,000 × 1,000 | — | 33.9 s (1 thread), 8.8 s (8 threads) |

**Memory (C++):** 39 MB on the plant data, where the Java analysis uses ~5 GB; 112 MB at
1,000 × 1,000, where Java needs ~6.7 GB (Supp. Table 5).

**Thread scaling on the plant data** (ms per sweep, dynamic chunks): 9.5 (1 thread), 5.7 (2),
4.6 (3), 4.1 (4), 3.1 (6), 2.8 (8). Splitting trees statically per thread regressed to 6.0 ms at
8 threads because of the slower efficiency cores and uneven polytomy cost.

**Other observations**

- About 40% of the time of an FM move on the plant data is spent outside the sweep (the taxon move
  and combining per-chunk results): 13.5 ms per move against 8 ms for the sweep.
- Each top-level FM pass keeps only about half of its moves: 96/200, 559/1,178 and 485/1,000.

## 5. Open items

- **ASTRAL-II estimated gene trees** (Dryad doi:10.6076/D10C7C, `estimated-gene-trees.tar.bz`):
  repeat sections 2 and 4 on the real conditions with `tools/phase0/run_phase0.sh` and
  `tools/compare/compare.sh`.
- **1KP compression numbers** from the run in progress.
- **Java at 1,000 × 1,000** on the same machine, and thread scaling on a Linux server with uniform
  cores.
