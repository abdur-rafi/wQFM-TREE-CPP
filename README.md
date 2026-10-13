# wqfm-tree-cpp

A C++20 reimplementation of [wQFM-TREE](https://github.com/abdur-rafi/wQFM-TREE) (Rafi, Rumi et al.,
*Bioinformatics Advances* 5(1), vbaf053, 2025), aiming to be much faster and lighter on memory, with
multithreading, while running the same search as the published Java version.

**Status: Phase 2 (complete search, optimized, single-threaded).** `build/wqfm-tree` runs the full
algorithm and produces the same species trees as the Java version on every test, including the
published 1KP plant tree. It is 42–110× faster than Java with about 1/50 of the memory
([docs/phase2-status.md](docs/phase2-status.md)). Parallelism and the native consensus come in Phase 3.
The repository also contains:

- `prototype/gainsweep.cpp`: the core of one FM iteration at the top level (score, gains of every
  taxon, taxon moves, one full FM pass). It reproduces the Java results **bit for bit** and runs
  19–78× faster single-threaded.
- Validation tools: brute-force quartet-enumeration checkers, and a harness that runs the Java and
  C++ code on the same input and compares every gain, every FM move and the final partition.
- Measurement tools: instrumentation for a copy of the Java code, a coalescent simulator, and a
  log analyzer.

See [docs/plan.md](docs/plan.md) for the design and phases, [docs/phase0-results.md](docs/phase0-results.md)
for the Phase 0 measurements, [docs/phase1-status.md](docs/phase1-status.md) for the Phase 1 validation, and
[docs/phase2-status.md](docs/phase2-status.md) for the Phase 2 optimizations.

## Layout

| Path | Contents |
|---|---|
| `src/` | the C++ implementation (`wqfm-tree`): parser, forest, scorer, consensus split, search, trace/replay |
| `tests/brute_check.cpp` | C++ scorer vs explicit quartet enumeration |
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
make            # build/wqfm-tree, build/brute_check, build/gainsweep
make java       # Java reference + instrumented builds (from ../src; override with WQFM_TREE_DIR=...)
make check      # brute-force checks + exact Java/C++ agreement; fails on any mismatch

# species tree (consensus tree as for the Java version, e.g. PAUP's greedy consensus)
build/wqfm-tree -i <gene trees> -c <consensus tree> -o <output> [--taxa-order <file>] [--trace <file>] [--no-compress]

# full-search validation against the Java version on your data
tools/compare/validate.sh <gene trees> <consensus tree> <out dir>

# Phase 0 measurement: compression, initial-split candidates, time split
tools/phase0/run_phase0.sh <gene trees> <out dir> [consensus tree]
```

`--taxa-order` takes `id label` lines or a Java trace. Taxon ids decide ties, and Java numbers taxa in
`HashSet` order, so pass Java's order to reproduce a Java run exactly (see the phase 1 notes on ties).

## Notes

- **exFAT drives.** macOS writes `._*` AppleDouble files next to every file it writes on exFAT.
  The scripts filter them out, and `.gitignore` ignores them. If git reports "non-monotonic index",
  delete `.git/objects/pack/._*`.
- **Upstream build quirks.** The upstream Java tree does not compile as-is (the legacy v1 classes
  are out of date). `tools/build_java.sh` excludes them; the algorithm (`Main` → `QFM2` → `*v2`
  classes) is unaffected.
