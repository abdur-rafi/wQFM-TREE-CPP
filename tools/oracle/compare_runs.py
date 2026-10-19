"""Compare two oracle sweeps (run_oracle.py), e.g. on estimated and on true gene trees.

usage: python3 -I tools/oracle/compare_runs.py ESTIMATED/results.csv TRUE/results.csv [--labels est,true]

For every condition and gene count present in both files (paired by replicate), prints:
  FN default      mean FN rate of the default runs on each set of gene trees, and the difference
  FN oracle       the same for the oracle runs (initial bipartitions from the true species tree)
  oracle effect   on each set: replicates where the oracle run's FN rate is lower / equal / higher
  >= true QS      replicates where the default estimate's quartet score is at least the true tree's
  FM left top     oracle runs in which FM left the best true-tree split at the top level
followed by the same quantities pooled by ILS level (species tree height) and over all conditions.
"""

import argparse
import csv
import statistics


def load(path):
    with open(path, newline='') as f:
        return {(r['condition'], r['replicate'], int(r['genes'])): r for r in csv.DictReader(f)}


def ils(cond):
    height = cond.split('.')[2]
    return {'500000': 'high ILS (500K)', '2000000': 'moderate ILS (2M)', '10000000': 'low ILS (10M)'}.get(height, height)


def stats(pairs):
    """pairs: list of (row on the first set, row on the second set)."""
    out = {'n': len(pairs)}
    for side in (0, 1):
        rs = [p[side] for p in pairs]
        fc = [float(r['fn_consensus']) for r in rs]
        fo = [float(r['fn_oracle']) for r in rs]
        out[side] = {
            'fn_cons': statistics.mean(fc),
            'fn_orcl': statistics.mean(fo),
            'orcl_better': sum(o < c for c, o in zip(fc, fo)),
            'orcl_worse': sum(o > c for c, o in zip(fc, fo)),
            'beats_true': sum(int(r['agree_consensus']) >= int(r['agree_true']) for r in rs),
            'top_left': sum(int(r['top_changed_oracle']) for r in rs),
        }
    return out


def line(name, st):
    n, a, b = st['n'], st[0], st[1]
    eff = lambda s: f'{s["orcl_better"]:3d}/{n - s["orcl_better"] - s["orcl_worse"]:3d}/{s["orcl_worse"]:3d}'
    return (f'{name:36s} {n:4d} | {a["fn_cons"]:.4f} {b["fn_cons"]:.4f} {b["fn_cons"] - a["fn_cons"]:+.4f} | '
            f'{a["fn_orcl"]:.4f} {b["fn_orcl"]:.4f} | {eff(a)}  {eff(b)} | '
            f'{100 * a["beats_true"] / n:4.0f}% {100 * b["beats_true"] / n:4.0f}% | '
            f'{100 * a["top_left"] / n:4.1f}% {100 * b["top_left"] / n:4.1f}%')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('first')
    ap.add_argument('second')
    ap.add_argument('--labels', default='est,true')
    a = ap.parse_args()
    l1, l2 = a.labels.split(',')
    r1, r2 = load(a.first), load(a.second)
    keys = sorted(set(r1) & set(r2))
    if not keys:
        raise SystemExit('no replicates in common')
    groups, pooled = {}, {}
    for k in keys:
        groups.setdefault((k[0], k[2]), []).append((r1[k], r2[k]))
        pooled.setdefault(ils(k[0]), []).append((r1[k], r2[k]))
        pooled.setdefault('all', []).append((r1[k], r2[k]))
    head = (f'{"condition, genes":36s} {"runs":>4s} | {"FN default":^22s} | {"FN oracle":^13s} | '
            f'{"oracle better/=/worse":^25s} | {">= true QS":^11s} | {"FM left top":^11s}')
    sub = (f'{"":36s} {"":4s} | {l1:>6s} {l2:>6s} {"diff":>7s} | {l1:>6s} {l2:>6s} | {l1:^12s} {l2:^12s} | '
           f'{l1:>5s} {l2:>5s} | {l1:>5s} {l2:>5s}')
    print(head)
    print(sub)
    print('-' * len(head))
    order = lambda g: (int(g[0].split('.')[1]), g[0], g[1])
    for (cond, genes) in sorted(groups, key=order):
        print(line(f'{cond}, {genes}', stats(groups[(cond, genes)])))
    print('-' * len(head))
    for name in ('low ILS (10M)', 'moderate ILS (2M)', 'high ILS (500K)', 'all'):
        if name in pooled:
            print(line(name, stats(pooled[name])))
    print(f'\nFN default/oracle: mean FN rate of the default and oracle runs on {l1} and {l2} gene trees '
          f'(diff = {l2} - {l1}).\noracle better/=/worse: replicates where the oracle run has a lower / equal / '
          f'higher FN rate than the default run.\n>= true QS: default estimates whose quartet score is at least '
          f'the true species tree\'s.\nFM left top: oracle runs where FM moved away from the best true-tree split '
          f'at the top level.')


if __name__ == '__main__':
    main()
