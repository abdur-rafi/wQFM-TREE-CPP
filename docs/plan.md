# Design and plan

## Principles

- **Same search as published.** Scoring, the consensus-based initial split, FM passes and the
  divide-and-conquer stay as in the paper and the Java code. Optimizations that change results
  (for example early termination of FM passes) are opt-in flags, evaluated separately.
- **Deterministic.** Identical output for any thread count: work is split into fixed chunks whose
  partial results are combined in a fixed order.
- **Verified at every step.** Brute-force quartet enumeration is the reference for scores and
  gains; a step-by-step comparison against the Java code covers the search itself.

## Where the Java time and memory go

- **Every FM move does two full passes over all gene trees.** Gains are recomputed for every taxon
  in O(N·d), and the swap then updates the cached "top branch" of every node in O(N). N is the
  total node count over all gene trees; d is the number of dummy taxa.
- **The initial split costs about as much as FM.** Each consensus candidate gets a full rebuild
  plus rescoring. It accounts for ~55% of measured runtime.
- **Every subproblem works on the full gene trees,** even with only a handful of taxa.
- **Polytomy nodes are expensive.** For a node of degree D, the Java code needs O(D²·d) to set up
  each subproblem, O(D²) per sweep for real-taxon gains, and O(d·D²) for dummy gains.
- **Memory.** Node keys are strings as long as the taxon list, kept for the whole run; every branch
  of every node holds a d-length array; JVM object overhead adds more.

## Design

1. **Flat storage.** All gene trees in one post-order array of 32-bit integers (parent, children
   offsets, taxon). Leaf positions per (tree, taxon) let a taxon move touch only its ancestors.
2. **Units.** Real and dummy taxa are both weighted "units", so one formula covers both:
   w(PA_ij) = a_i a_j − Σ_U u_i u_j and w(PB_k) = (b_k² − Σ_U u_k²)/2.
   Nodes store only below-subtree totals; the top branch is derived from per-tree totals.
3. **Compressed trees per subproblem.** Collapse every maximal subtree whose taxa all lie under one
   dummy taxon into a single weighted leaf, and merge same-dummy sibling leaves. This is exact,
   since such nodes only anchor quartets that the score excludes; `CompressCheck` verifies it.
   Each child subproblem's trees are derived from its parent's.
4. **Closed-form polytomy formulas.** Power sums give w(S(g,u)), w(U(g,u)) and every real-taxon
   gain of a degree-D node in O(D·(1+d)²), or O(D) with no dummy taxa: each dummy's correction is
   a rank-1 term, so f = w(S) + w(U)/2 is a function of additive per-branch aggregates.
   Implemented in Phase 1 (`PolyEval`), checked by brute force, and bit-exact against Java at the
   top level of the 1KP data (polytomies of up to 576 branches).
5. **Dummy-taxon gains.** A low-degree polynomial in the dummy's per-tree weight, plus a sparse
   correction on the nodes above that dummy's leaves.
6. **Node and tree hashing.** 128-bit hashes for duplicate detection, including collapsed-leaf
   weights; no long string keys. Within a subproblem, identical compressed gene trees can be
   stored once with a count.
7. **Integer arithmetic at the top level.** With no dummy taxa every quantity is an integer, so
   results there are exact.
8. **Native greedy consensus and Newick parser.** Removes the PAUP, Perl and Python
   dependencies; a given consensus file can still be supplied for exact reproduction.
9. **Parallelism.** Within a subproblem, gene trees are processed in fixed chunks claimed
   dynamically by threads (static splitting regressed on heterogeneous cores). Subproblems run
   as parallel tasks only when large enough to be worth it, with a cap on live tasks to bound
   memory. After compression, the top levels dominate, so parallelism within a subproblem is
   what matters most.
10. **Numerics.** Doubles; no `-ffast-math`; `-ffp-contract=off`; Java's EPS = 1e-5 rule for
    cumulative gains.

## Phases

| Phase | Work |
|---|---|
| 0 | Measurements on simulated and real data; brute-force checkers; top-level prototype. **Done, except real estimated gene trees.** |
| 1 | Single-threaded C++ core: flat storage, unit formulation, full divide-and-conquer with dummy taxa, consensus initial split, closed-form polytomies, plain O(d) dummy loop. Validated against the checkers and Java traces. **Implemented; same trees as Java on all simulated tests ([phase1-status.md](phase1-status.md)); 1KP comparison in progress.** |
| 2 | Optimizations, each verified: compression, sparse dummy gains, hashing, integer top level, faster sweep kernel. **Compression (derived from the parent's forest), single setup, incremental candidates and the no-dummy polytomy fast path are done: 4–7× faster than Phase 1, 42–110× faster than Java, same trees ([phase2-status.md](phase2-status.md)). Sparse dummy gains and integer top level measured as not worth it; duplicate-tree storage and candidate batching deferred.** |
| 3 | Native consensus and Newick parser (fuzz-tested); parallelism at both levels with deterministic results. |
| 4 | Benchmarks against Java, ASTRAL and TREE-QMC (ASTRAL-II conditions, 2000-taxon, 1KP, avian); CMake, static binary, bioconda; release. |

**Optional later:** early termination of FM passes (each top-level pass keeps only about half of
its moves); quartet-based branch support from the same machinery.

## Corrections made during planning

- **The FM step is inherently Θ(N) per move,** because every taxon's gain changes at every node.
  The redesign removes the d factor and the constant overhead, not N itself.
- **Deduplicating candidate splits is worth only ~10%.** Java's size filter already limits each
  subproblem to about 2s candidates (s = subproblem size), so the paper's O(n³k log n) bound is
  loose; the actual behaviour is about O(n²k log n).
