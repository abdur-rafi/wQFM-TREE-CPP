"""Compare the FN rates of two sweeps, paired by (condition, replicate, genes).

usage: python3 -I tools/weighted/compare_fn.py BASE.csv OTHER.csv [--labels base,other] [--col-a fn_rate] [--col-b fn_rate]

Accepts run_all.py results (column fn_rate) and run_oracle.py results (fn_consensus / fn_oracle);
pass --col-a / --col-b to pick the column of each file. Prints, per condition and gene count and
pooled by ILS level (species tree height), the mean FN of both sweeps, the change, how many
replicates are better / equal / worse with OTHER, and the two-sided sign-test p-value.
"""

import argparse
import csv
import math
import statistics


def load(path, col):
    out = {}
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            if r.get(col, '') != '':
                out[(r['condition'], r['replicate'], int(r['genes']))] = float(r[col])
    return out


def sign_p(w, l):
    n = w + l
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(min(w, l) + 1)) / 2 ** n)


def ils(cond):
    return {'500000': 'high ILS', '2000000': 'moderate ILS', '10000000': 'low ILS'}.get(cond.split('.')[2], cond)


def line(name, pairs):
    a = [x for x, _ in pairs]
    b = [y for _, y in pairs]
    w = sum(y < x for x, y in pairs)
    l = sum(y > x for x, y in pairs)
    ma, mb = statistics.mean(a), statistics.mean(b)
    return (f'{name:38s} {len(pairs):4d} | {ma:.4f} {mb:.4f} {mb - ma:+.4f} ({100 * (mb - ma) / ma if ma else 0:+5.1f}%) | '
            f'{w:4d} {len(pairs) - w - l:4d} {l:4d} {sign_p(w, l):8.1e}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('base')
    ap.add_argument('other')
    ap.add_argument('--labels', default='base,other')
    ap.add_argument('--col-a', default=None)
    ap.add_argument('--col-b', default=None)
    a = ap.parse_args()
    la, lb = a.labels.split(',')
    guess = lambda p: 'fn_consensus' if 'oracle' in p else 'fn_rate'
    A, B = load(a.base, a.col_a or guess(a.base)), load(a.other, a.col_b or guess(a.other))
    keys = sorted(set(A) & set(B), key=lambda k: (int(k[0].split('.')[1]), k[0], k[2], k[1]))
    if not keys:
        raise SystemExit('no replicates in common')
    head = f'{"condition, genes":38s} {"n":>4s} | {la:>6s} {lb:>6s} {"change":>16s} | {"better":>4s} {"=":>4s} {"worse":>4s} {"p":>8s}'
    print(head)
    print('-' * len(head))
    groups, pooled = {}, {}
    for k in keys:
        groups.setdefault((k[0], k[2]), []).append((A[k], B[k]))
        pooled.setdefault(ils(k[0]), []).append((A[k], B[k]))
        pooled.setdefault(f'{k[2]} genes', []).append((A[k], B[k]))
        pooled.setdefault('all', []).append((A[k], B[k]))
    for (cond, genes), pairs in groups.items():
        print(line(f'{cond}, {genes}', pairs))
    print('-' * len(head))
    for name in ['low ILS', 'moderate ILS', 'high ILS'] + sorted((x for x in pooled if x.endswith('genes')), key=lambda x: int(x.split()[0])) + ['all']:
        if name in pooled:
            print(line(name, pooled[name]))
    print(f'\nbetter / = / worse: replicates where {lb} has a lower / equal / higher FN rate than {la}; p: two-sided sign test.')


if __name__ == '__main__':
    main()
