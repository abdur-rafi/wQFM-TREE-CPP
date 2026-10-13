# wqfm-tree-cpp

A C++20 reimplementation of [wQFM-TREE](https://github.com/abdur-rafi/wQFM-TREE) (Rafi, Rumi et al.,
*Bioinformatics Advances* 5(1), vbaf053, 2025), aiming to be much faster and lighter on memory, with
multithreading, while running the same search as the published Java version.

**Status: Phase 0 (measurement and prototype).** This is not yet a usable species-tree tool. It contains:

- `prototype/gainsweep.cpp`: the core of one FM iteration at the top level (score, gains of every
  taxon, taxon moves, one full FM pass). It reproduces the Java results **bit for bit** and runs
  19–78× faster single-threaded.
- Validation tools: brute-force quartet-enumeration checkers, and a harness that runs the Java and
  C++ code on the same input and compares every gain, every FM move and the final partition.
- Measurement tools: instrumentation for a copy of the Java code, a coalescent simulator, and a
  log analyzer.

See [docs/plan.md](docs/plan.md) for the design and phases, and
[docs/phase0-results.md](docs/phase0-results.md) for the measurements so far.

## Layout

| Path | Contents |
|---|---|
| `prototype/gainsweep.cpp` | top-level gain sweep and FM pass (flat arrays, O(D) polytomy formulas, thread pool) |
| `tools/java/` | `TopLevelBench` (Java side of the comparison), `BruteCheck`, `CompressCheck`, `Phase0Stats`, `instrument.py` |
| `tools/build_java.sh` | builds the Java reference and instrumented copies from the parent wQFM-TREE `src/` |
| `tools/check.sh` | all correctness checks (`make check`) |
| `tools/compare/` | Java vs C++ comparison on any gene-tree file |
| `tools/phase0/` | Phase 0 measurement on any gene-tree file, and the log analyzer |
| `tools/greedy_consensus.py` | PAUP-free greedy consensus (complete gene trees) and Newick cleaning |
| `tools/sim/` | multispecies-coalescent simulator used for the synthetic measurements |

This repository lives inside a wQFM-TREE checkout and is ignored by it. The Java tooling compiles
from a copy of `../src`, so the upstream sources are never modified.

## Build and run

Requirements: a C++20 compiler with `<barrier>` (Apple clang 16, GCC ≥ 11, clang ≥ 14), a JDK
(tested with 21), and Python 3.

```sh
make            # build/gainsweep
make java       # Java reference + instrumented builds (from ../src; override with WQFM_TREE_DIR=...)
make check      # brute-force checks + exact Java/C++ agreement; fails on any mismatch

# Java vs C++ on your own gene trees (one sweep + one full FM pass, same partition)
tools/compare/compare.sh <gene trees> <out dir> [threads]

# Phase 0 measurement: compression, initial-split candidates, time split
tools/phase0/run_phase0.sh <gene trees> <out dir> [consensus tree]
```

The prototype CLI:
`build/gainsweep <gene trees> <taxa order> <partition> <threads> <reps> sweep|fm <out prefix>`.
The taxa-order and partition files come from `TopLevelBench`, so taxon ids and tie-breaking match
the Java version.

## Notes

- **exFAT drives.** macOS writes `._*` AppleDouble files next to every file it writes on exFAT.
  The scripts filter them out, and `.gitignore` ignores them. If git reports "non-monotonic index",
  delete `.git/objects/pack/._*`.
- **Upstream build quirks.** The upstream Java tree does not compile as-is (the legacy v1 classes
  are out of date). `tools/build_java.sh` excludes them; the algorithm (`Main` → `QFM2` → `*v2`
  classes) is unaffected.
