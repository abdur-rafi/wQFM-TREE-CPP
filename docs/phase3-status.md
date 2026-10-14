# Phase 3 status

Multithreading and a built-in greedy consensus.

## Threads (`-t N`; the default was all hardware threads, now 1)

- **Within a subproblem.** The gene trees of each subproblem's forest are split into chunks of
  about 4,096 nodes (at most 128). Their number depends on the forest only. Score, gain sweeps,
  partition setup and dummy-taxon moves run over the chunks in parallel, and per-chunk results are
  summed in chunk order.
- **Across subproblems.** After each divide step, child 1 runs as a task picked up by an idle
  thread and child 0 runs on the current thread (fork-join).
- **The thread pool** (`src/thread_pool.*`) supports nested, concurrent `parallel_for` calls and
  tasks from any thread. Items are claimed dynamically, which matters on CPUs with fast and slow
  cores. A thread waiting for a task runs queued tasks and loop items meanwhile.
- **Determinism.** Each subproblem is computed the same way whatever thread runs it, so the output
  is byte-identical for any thread count. Dummy-taxon ids depend on scheduling, but they only label
  leaves and never affect a computation. Tracing and replay process subproblems one at a time
  (their per-tree work stays parallel).
- **Shared state.** The dummy-taxon pool has reserved storage and a lock. The spare-scorer list and
  timers are locked; counters are atomic.

**Checks.** ThreadSanitizer and AddressSanitizer runtimes crash on this machine (macOS 26.4; even
trivial programs fail), and GCC 15's headers do not build against the SDK. Instead:
- every piece of shared state was reviewed (pool job/task lifetimes, publication of new dummy
  taxa, per-thread scratch, disjoint chunk writes);
- 90 runs at 3, 5 and 8 threads on three datasets were all byte-identical to the 1-thread output;
- replay of Java traces with 8 threads agrees with Java;
- `make check` compares 1 and 4 threads.

**To do:** run ThreadSanitizer on Linux.

The review also fixed an error-path issue: if a child subproblem throws, the parent now waits for
the sibling task, which refers to its stack frame, before unwinding.

## Built-in greedy consensus (`-c` now optional)

`src/greedy.*`:
- **Counting splits:** gene-tree splits are counted with 128-bit XOR hashes, and the two edges
  below a degree-2 root are counted as one.
- **Building the tree:** splits are added in order of decreasing frequency (ties: first seen),
  keeping each split compatible with those already accepted. Compatibility is tested at the lowest
  common ancestor in the growing tree, and the loop stops once the tree is fully resolved.
- **Complete gene trees only.** Otherwise it fails with a message asking for `-c`, since PAUP's
  handling of incomplete trees is unknown.
- `--write-consensus <file>` saves the consensus used.

Results:
- **11-taxon example:** identical to PAUP's greedy consensus (RF = 0).
- **Simulated sets:** identical to the Python reference implementation on two of three. On the
  300-taxon set, 13 splits differ through frequency ties broken in a different order.
- **Species trees.** Search results with the built-in consensus match Java (run on the reference
  consensus) on two of three simulated sets. On the third the consensus has the same splits, but
  its Newick layout changes the order in which initial-split candidates are scored, so exactly-tied
  candidates break differently (2 splits differ). Pass the PAUP consensus with `-c` to reproduce a
  Java run exactly.

## Timings

Apple M3 (4 performance + 4 efficiency cores), idle machine, minimum of 3 runs (1KP: single runs).
All runs give identical trees for every thread count.

| Dataset | Taxa × genes | Java | 1 thread | 4 threads | 8 threads | 8 threads vs Java |
|---|---|---|---|---|---|---|
| Simulated, discordance 0.44 | 200 × 200 | 44.6 s | 0.86 s | 0.35 s | 0.33 s | 135× |
| Simulated, discordance 0.77 | 200 × 200 | 82.6 s | 1.15 s | 0.45 s | 0.41 s | 200× |
| 15% missing taxa + polytomies | 120 × 150 | 25.9 s | 0.40 s | 0.17 s | 0.19 s | 136× |
| Simulated, discordance 0.34 | 300 × 150 | 84.6 s | 1.41 s | 0.60 s | 0.57 s | 148× |
| 1KP plant | 1,178 × 410 | ≈5.5 h (paper) | 153 s | – | 42.4 s | ≈470× |

Java times for the simulated sets come from the traced validation runs on the same machine. The
1KP output is the published wQFM-TREE tree at every thread count. Peak memory on the 1KP data is
222 MB with 1 thread and 260 MB with 8 (more scorers alive at once); the Java version needs ~5.7 GB.

Going from 4 to 8 threads gains little here because the extra 4 cores are efficiency cores; a
server with uniform cores should scale further. Not yet measured.
