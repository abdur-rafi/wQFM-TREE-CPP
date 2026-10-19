# Oracle experiment: initial bipartitions from the true species tree

**Question.** Where should algorithmic work go: the initial bipartition, or the scoring (the
objective FM optimises)? If wQFM-TREE's estimates are limited by poor starting points, starting
every subproblem from the true species tree's bipartitions should give better trees.

## Method

Every replicate is run twice on the same gene trees, with the same build and options:

| Run | Initial bipartition of each subproblem |
|---|---|
| `consensus` (default) | the best-scoring edge of the greedy consensus of the gene trees |
| `oracle` | the best-scoring edge of the **true species tree** (`-c <true tree>`) |

In both runs, the FM passes, divide step and dummy taxa are unchanged. Edges are projected onto a
subproblem the usual way: a real taxon by its side, a dummy taxon to the side holding at least half
its weight. At the top level (no dummy taxa) the oracle run therefore starts from the true-tree edge
with the highest wQFM score.

Each estimate, and the true tree, is scored by:

- **Quartet score.** The fraction of gene-tree quartets the tree agrees with: the score ASTRAL
  reports, computed exactly by `build/qscore` and checked against explicit enumeration in
  `brute_check`.
- **FN rate.** The fraction of true-tree branches missing from the estimate.

Each run also writes a subproblem log (`wqfm-tree --subproblem-log`). For each subproblem it records
the score of the initial bipartition, the score after FM, and whether FM changed the bipartition. In
the oracle run, a change at the top level means the wQFM objective scores a bipartition that is not in
the true tree above every true-tree edge.

Reading the results:

- **Oracle better than consensus** (higher quartet score or lower FN rate): better initial
  bipartitions would help.
- **Oracle and consensus reach the same trees:** FM is not limited by its starting point.
- **Default estimates score at least as high as the true tree:** maximising the quartet score
  harder cannot make them more accurate. Gains must then come from the objective.
- **FM moves the oracle off true-tree edges:** the wQFM objective itself prefers wrong splits.

## Reproducing

Data: the ASTRAL-II simulated datasets (Mirarab and Warnow 2015) in the `run_all.py` layout:

```
DATA/estimated-gene-trees/<condition>/<replicate>/estimatedgenetre
DATA/true_species_trees/<condition>/<replicate>/s_tree.trees
```

```sh
make                                   # build/wqfm-tree and build/qscore
build/brute_check 1 300                # includes the quartet-score check
python3 -I tools/oracle/run_oracle.py -i DATA -o OUT -t 8 --genes 50,200,1000
```

The nine conditions with at least 100 taxa run by default (`--conditions` selects others). The
threads only change the speed. The trees, scores and subproblem logs, apart from the line order of
the logs, are the same for any `-t`. Runs are resumable: rerunning the command skips finished jobs.
`--overwrite` reruns them.

Outputs in `OUT`:

| File | Contents |
|---|---|
| `results.csv` | one row per replicate and gene count: quartet scores (true, consensus, oracle), exact quartet counts, FN/FP rates, RF distance between the two estimates, FM statistics, run times |
| `summary.txt` | per condition and gene count: means, win/tie/loss counts with sign-test p-values, FM statistics |
| `setup.json` | the command, repository commit, SHA-256 of both executables, and the date of each invocation |
| `<condition>/<replicate>/<genes>-genes/` | `consensus/` and `oracle/` (species tree, log, subproblem log), `qscore.tsv` |

`python3 -I tools/oracle/run_oracle.py -o OUT --summary-only` rebuilds `summary.txt` from
`results.csv`.

## Results

Pending: the full sweep (9 conditions × 50 replicates × 50, 200 and 1000 genes = 1,350 jobs of two
runs each) is running.
