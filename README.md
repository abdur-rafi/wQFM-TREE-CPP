# wqfm-tree-cpp

A C++20 reimplementation of [wQFM-TREE](https://github.com/abdur-rafi/wQFM-TREE) (Rafi, Rumi et al.,
*Bioinformatics Advances* 5(1), vbaf053, 2025), aiming to be much faster and lighter on memory, with
multithreading, while running the same search as the published Java version.

**Status: Phase 3 done; ASTRAL-II benchmark done.**
- `build/wqfm-tree` runs the full algorithm without PAUP, and its output is byte-identical for any
  thread count.
- **Same trees as Java.** It matches the Java version on every test, including the published 1KP
  plant tree. On the paper's ASTRAL-II datasets, Java-trace replays agree number by number, and over
  1,101 runs accuracy equals the published trees.
- **Speed and memory.** On an 8-core laptop it is 110–240× faster than the Java times in the paper,
  with about 1/6–1/9 of the memory.

See [docs/astral2-results.md](docs/astral2-results.md) and [docs/phase3-status.md](docs/phase3-status.md).
The repository also contains:

- `prototype/gainsweep.cpp`: the core of one FM iteration at the top level (score, gains of every
  taxon, taxon moves, one full FM pass). It reproduces the Java results **bit for bit** and runs
  19–78× faster single-threaded.
- Validation tools: brute-force quartet-enumeration checkers, and a harness that runs the Java and
  C++ code on the same input and compares every gain, every FM move and the final partition.
- Measurement tools: instrumentation for a copy of the Java code, a coalescent simulator, and a
  log analyzer.

See [docs/plan.md](docs/plan.md) for the design and phases, [docs/phase0-results.md](docs/phase0-results.md)
for the Phase 0 measurements, [docs/phase1-status.md](docs/phase1-status.md) for the Phase 1 validation,
[docs/phase2-status.md](docs/phase2-status.md) for the Phase 2 optimizations,
[docs/phase3-status.md](docs/phase3-status.md) for multithreading and the built-in consensus, and
[docs/astral2-results.md](docs/astral2-results.md) for validation, accuracy and speed on the paper's
ASTRAL-II datasets (and the kernel work done there), and [docs/oracle-experiment.md](docs/oracle-experiment.md)
for the oracle experiment (initial bipartitions from the true species tree).

## Layout

| Path | Contents |
|---|---|
| `src/` | the C++ implementation (`wqfm-tree`): parser, forest, scorer, consensus split, search, thread pool, greedy consensus, trace/replay; `qscore` (exact quartet score of a species tree) |
| `tests/brute_check.cpp` | C++ scorer, polytomy resolution and quartet score vs explicit quartet enumeration |
| `tests/bench_scorer.cpp` | kernel micro-benchmark (ns per node, checksum of all results) |
| `tools/astral2/` | ASTRAL-II runs scored against true and published trees; Java/C++ validation batches; `run_methods.py`: weighted ASTRAL, TREE-QMC or wQFM-TREE sweeps (parallel single-thread jobs, resumable) |
| `tools/oracle/` | oracle experiment: default vs true-tree initial bipartitions, quartet scores and FN rates |
| `tools/weighted/` | low-support contraction baseline; `compare_fn.py`: paired FN comparison of two sweeps (sign test) |
| `prototype/gainsweep.cpp` | Phase 0 top-level gain sweep and FM pass (multithreaded) |
| `tools/java/` | `TopLevelBench`, `BruteCheck`, `CompressCheck`, `Phase0Stats` + `instrument.py`, `Trace` + `trace_instrument.py` |
| `tools/build_java.sh` | builds the Java reference and instrumented copies from the parent wQFM-TREE `src/` |
| `tools/check.sh` | all correctness checks (`make check`) |
| `tools/compare/` | `validate.sh` (full search: Java trace vs C++ run and replay), `trace_compare.py`, `rf.py`, `compare.sh` (prototype) |
| `tools/phase0/` | Phase 0 measurement on any gene-tree file, and the log analyzer |
| `tools/greedy_consensus.py` | PAUP-free greedy consensus (complete gene trees) and Newick cleaning |
| `tools/sim/` | multispecies-coalescent simulator used for the synthetic measurements |

This repository lives inside a wQFM-TREE checkout and is ignored by it. The Java tooling compiles
from a copy of `../src`, so the upstream sources are never modified.

## Build and run

Requirements: a C++20 compiler with `<barrier>` (Apple clang 16, GCC ≥ 11, clang ≥ 14), a JDK
(tested with 21), and Python 3.

```sh
make            # build/wqfm-tree, build/qscore, build/brute_check, build/bench_scorer, build/gainsweep
make java       # Java reference + instrumented builds (from ../src; override with WQFM_TREE_DIR=...)
make check      # brute-force checks + exact Java/C++ agreement; fails on any mismatch

# species tree (built-in greedy consensus; or pass one with -c, e.g. PAUP's as in run.sh)
build/wqfm-tree -i <gene trees> -o <output> [-c <consensus tree>] [-t <threads>]
                [--write-consensus <file>] [--relative-eps] [--resolve-polytomies] [--taxa-order <file>]
                [--trace <file>] [--no-compress] [--subproblem-log <file>]

# quartet score of species trees (exact; ASTRAL's normalised score is the third column)
build/qscore -g <gene trees> [-t <threads>] <species tree file>...

# full-search validation against the Java version on your data
tools/compare/validate.sh <gene trees> <consensus tree> <out dir>

# Phase 0 measurement: compression, initial-split candidates, time split
tools/phase0/run_phase0.sh <gene trees> <out dir> [consensus tree]
```

`--taxa-order` takes `id label` lines or a Java trace. Taxon ids decide ties, and Java numbers taxa in
`HashSet` order, so pass Java's order to reproduce a Java run exactly (see the phase 1 notes on ties).

## Notes

- **Threads.** `-t` defaults to 1. Every FM move ends with a parallel loop, so when several jobs
  share a machine, keep their total threads within the core count: oversubscribed runs slow down
  ~10×, not 2×.
- **`--relative-eps`.** Java keeps a longer prefix of FM moves only if its cumulative gain beats the
  best by more than an absolute 1e-5. At large score scales, rounding in the running sum can exceed
  that (see [docs/astral2-results.md](docs/astral2-results.md)); this option uses
  1e-5 + 1e-12 × (sum of |gain| so far). Off by default, so results match the Java version.
- **`--resolve-polytomies`.** Resolves gene-tree polytomies arbitrarily before the search, replacing
  the upstream `runResolvingPolytomy.sh` step. It uses DendroPy's deterministic rule
  (`resolve_polytomies()` without a random generator): join a node's first two children under a new
  node, appended as the last child, until two children remain, so the resolved trees are exactly
  DendroPy's. As in that script, the built-in consensus is still built from the unresolved trees.
  Polytomies are scored exactly without this option, so it is not needed for speed.

- **exFAT drives.** macOS writes `._*` AppleDouble files next to every file it writes on exFAT.
  The scripts filter them out, and `.gitignore` ignores them. If git reports "non-monotonic index",
  delete `.git/objects/pack/._*`.
- **Upstream build quirks.** The upstream Java tree does not compile as-is (the legacy v1 classes
  are out of date). `tools/build_java.sh` excludes them; the algorithm (`Main` → `QFM2` → `*v2`
  classes) is unaffected.
