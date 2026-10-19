"""Compare variants of the wQFM-TREE objective (wqfm-tree --kappa / --power / --dummy-gamma).

Every variant runs on every selected replicate from up to three starts (--arms): the default initial
bipartitions (greedy consensus), the true species tree's ("oracle", -c <true tree>), and the tree of
the published objective's default run ("seeded", -c OUT/base/default/.../species.tre; a two-stage
search). From the oracle start, a better objective keeps more true splits; the default start checks
that a variant is not merely reluctant to move. Results are paired by replicate against the published
objective ("base") from the same start, and seeded runs against base's default run.

usage: python3 -I tools/oracle/objectives.py -i DATA -o OUT [-j JOBS] [-t THREADS] [options]

  -i DATA               data in the run_all.py layout (true species trees in DATA/true_species_trees)
  -o OUT                results: OUT/results.csv, OUT/summary.txt, OUT/<variant>/<arm>/<cond>/<rep>/<genes>/
  --gene-trees-dir, --gene-file   gene trees (default: true-gene-trees, truegenetrees)
  --conditions LIST     default: the six 200-taxon conditions
  --replicates LIST     default: 1-20
  --genes LIST          default: 50,200,1000
  --variants LIST       names from VARIANTS (default: all)
  --arms LIST           starts: default,oracle,seeded (default: default,oracle)
  -j JOBS               runs at a time (default 4); -t THREADS per run (default 2)
  --summary-only        rebuild summary.txt from results.csv

Variants: score per subproblem = (2F - kappa T) / T^power, F = satisfied + unresolved/2, T = all
quartets with two units on each side; dummy taxa of m members weigh m^gamma in total.
  base               kappa 1, power 0, gamma 0   satisfied - violated (published)
  alpha<a>           kappa 2a/(1+a)              satisfied - a * violated
  power<p>           power p                     (satisfied - violated) / T^p; p = 0.5: a z-score,
                                                 p = 1: the fraction satisfied, rescaled
  gamma<g>           gamma g                     dummy taxa weigh (members)^g instead of 1
  init-alpha0.5      alpha 1/2 for choosing the initial bipartition only, FM with alpha 1
  fm-alpha0.5        alpha 1/2 for the FM passes only, initial bipartition chosen with alpha 1
  wsupport, wlength, whybrid   weighted quartets (--weight support | length | hybrid)
"""

import argparse
import concurrent.futures as cf
import csv
import math
import os
import statistics
import subprocess
import sys
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HERE, '..', 'astral2'))
sys.path.insert(0, os.path.join(HERE, '..', 'compare'))
import rf  # noqa: E402
import run_all  # noqa: E402

sys.setrecursionlimit(100000)


def alpha(a):
    return {'kappa': 2 * a / (1 + a)}


VARIANTS = {
    'base': {},
    'alpha0.5': alpha(0.5), 'alpha0.75': alpha(0.75), 'alpha1.5': alpha(1.5), 'alpha2': alpha(2), 'alpha3': alpha(3),
    'power0.25': {'power': 0.25}, 'power0.5': {'power': 0.5}, 'power1': {'power': 1},
    'gamma0.5': {'dummy_gamma': 0.5}, 'gamma1': {'dummy_gamma': 1},
    # alpha 1/2 in one place only: choosing the initial bipartition, or the FM passes
    'init-alpha0.5': {'init_kappa': 2 / 3}, 'fm-alpha0.5': {'kappa': 2 / 3, 'init_kappa': 1.0},
    # weighted quartets (wqfm-tree --weight), as weighted ASTRAL's support, length and hybrid weighting
    'wsupport': {'weight': 'support'}, 'wlength': {'weight': 'length'}, 'whybrid': {'weight': 'hybrid'},
}
ARMS = ('default', 'oracle', 'seeded')
CONDITIONS = ['model.200.10000000.0.0000001', 'model.200.10000000.0.000001', 'model.200.2000000.0.0000001',
              'model.200.2000000.0.000001', 'model.200.500000.0.0000001', 'model.200.500000.0.000001']
FIELDS = ['variant', 'arm', 'condition', 'replicate', 'genes', 'fn', 'top_changed']


def flags(v):
    out = []
    for k, flag in (('kappa', '--kappa'), ('power', '--power'), ('dummy_gamma', '--dummy-gamma'),
                    ('init_kappa', '--init-kappa'), ('weight', '--weight')):
        if k in VARIANTS[v]:
            x = VARIANTS[v][k]
            out += [flag, x if isinstance(x, str) else repr(x)]
    return out


def top_changed(log):
    with open(log) as f:
        for line in f:
            if not line.startswith('#') and line.split('\t')[1] == '0':
                return int(line.split('\t')[5])
    return 0


def fn_rate(true_newick, est_newick):
    st, tt = rf.splits(rf.parse(true_newick))
    se, te = rf.splits(rf.parse(est_newick))
    if tt != te:
        raise RuntimeError('estimate and true tree have different taxa')
    return len(st - se) / len(st)


def sign_p(w, l):
    n = w + l
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(min(w, l) + 1)) / 2 ** n)


def ils(cond):
    return {'500000': 'high', '2000000': 'moderate', '10000000': 'low'}.get(cond.split('.')[2], '?')


def summary(rows, path):
    by = {(r['variant'], r['arm'], r['condition'], r['replicate'], r['genes']): r for r in rows}
    variants = [v for v in VARIANTS if any(k[0] == v for k in by)]
    lines = []
    for arm in [x for x in ARMS if any(k[1] == x for k in by)]:
        ref = 'default' if arm == 'seeded' else arm       # seeded runs are compared with base's default run
        lines.append(f'== {arm} start: mean FN rate, and FN vs base from the {ref} start '
                     f'(better/equal/worse, sign test p, mean change)')
        groups = [('all', lambda k: True)] + [(f'{x} ILS', lambda k, x=x: ils(k[2]) == x) for x in ('low', 'moderate', 'high')] \
            + [(f'{g} genes', lambda k, g=g: k[4] == g) for g in sorted({k[4] for k in by}, key=int)]
        head = f'{"variant":11s} ' + ' '.join(f'{g:>9s}' for g, _ in groups) + ' | vs base: better/=/worse      p   mean dFN'
        lines.append(head)
        for v in variants:
            keys = [k for k in by if k[0] == v and k[1] == arm]
            cells = []
            for _, sel in groups:
                fs = [float(by[k]['fn']) for k in keys if sel(k)]
                cells.append(f'{statistics.mean(fs):9.4f}' if fs else f'{"-":>9s}')
            refk = lambda k: ('base', ref) + k[2:]
            pairs = [(float(by[refk(k)]['fn']), float(by[k]['fn'])) for k in keys if refk(k) in by]
            w = sum(b > x for b, x in pairs)
            l = sum(b < x for b, x in pairs)
            d = statistics.mean(x - b for b, x in pairs) if pairs else 0.0
            lines.append(f'{v:11s} ' + ' '.join(cells) + f' | {w:4d}/{len(pairs) - w - l:4d}/{l:4d} {sign_p(w, l):9.1e} {d:+.5f}')
        top = [(v, statistics.mean(int(by[k]['top_changed']) for k in by if k[0] == v and k[1] == 'oracle'))
               for v in variants if any(k[0] == v and k[1] == 'oracle' for k in by)]
        lines.append('')
    lines.append('oracle start, top-level true split left by FM: ' + ', '.join(f'{v} {100 * t:.1f}%' for v, t in top))
    text = '\n'.join(lines) + '\n'
    run_all.write_atomic(path, text)
    print(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-i', '--input')
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('--gene-trees-dir', default='true-gene-trees')
    ap.add_argument('--gene-file', default='truegenetrees')
    ap.add_argument('--conditions', default=','.join(CONDITIONS))
    ap.add_argument('--replicates', default='1-20')
    ap.add_argument('--genes', default='50,200,1000')
    ap.add_argument('--variants', default=','.join(VARIANTS))
    ap.add_argument('--arms', default='default,oracle')
    ap.add_argument('-j', '--jobs', type=int, default=4)
    ap.add_argument('-t', '--threads', type=int, default=2)
    ap.add_argument('--binary', default=os.path.join(REPO, 'build', 'wqfm-tree'))
    ap.add_argument('--summary-only', action='store_true')
    a = ap.parse_args()
    out = os.path.abspath(a.output)
    csv_path = os.path.join(out, 'results.csv')
    if a.summary_only:
        with open(csv_path, newline='') as f:
            summary(list(csv.DictReader(f)), os.path.join(out, 'summary.txt'))
        return
    data = os.path.abspath(a.input)
    os.makedirs(out, exist_ok=True)
    variants = a.variants.split(',')
    for v in variants:
        if v not in VARIANTS:
            sys.exit(f'unknown variant {v}; known: {", ".join(VARIANTS)}')
    reps = [f'{r:02d}' for r in run_all.parse_list(a.replicates)]
    genes = run_all.parse_list(a.genes)
    rows = {}
    if os.path.exists(csv_path):
        with open(csv_path, newline='') as f:
            for r in csv.DictReader(f):
                rows[(r['variant'], r['arm'], r['condition'], r['replicate'], r['genes'])] = r
    jobs = []
    for cond in a.conditions.split(','):
        for rep in reps:
            src = os.path.join(data, a.gene_trees_dir, cond, rep, a.gene_file)
            true_path = os.path.join(data, run_all.TRUE_TREES_DIR, cond, rep, 's_tree.trees')
            if not os.path.isfile(src) or not os.path.isfile(true_path):
                print(f'skip {cond}/{rep}: missing input', file=sys.stderr)
                continue
            for g in genes:
                for v in variants:
                    for arm in a.arms.split(','):
                        if (v, arm, cond, rep, str(g)) not in rows:
                            jobs.append((v, arm, cond, rep, g, src, true_path))
    print(f'{len(jobs)} runs to do ({len(rows)} done); {a.jobs} at a time, {a.threads} threads each', flush=True)
    lock = threading.Lock()
    cut = {}                                    # the first g gene trees, written once per (cond, rep, g)

    def gene_file(cond, rep, g, src):
        key = (cond, rep, g)
        with lock:
            if key not in cut:
                p = os.path.join(out, 'genes', cond, rep, f'{g}.tre')
                if not os.path.isfile(p):
                    os.makedirs(os.path.dirname(p), exist_ok=True)
                    run_all.write_atomic(p, '\n'.join(run_all.read_trees(src)[:g]) + '\n')
                cut[key] = p
            return cut[key]

    def run(job):
        v, arm, cond, rep, g, src, true_path = job
        d = os.path.join(out, v, arm, cond, rep, f'{g}-genes')
        os.makedirs(d, exist_ok=True)
        cmd = [a.binary, '-i', gene_file(cond, rep, g, src), '-o', os.path.join(d, 'species.tre'), '-t', str(a.threads),
               '--quiet', '--subproblem-log', os.path.join(d, 'subproblems.tsv')] + flags(v)
        if arm == 'oracle':
            cmd += ['-c', true_path]
        elif arm == 'seeded':
            seed = os.path.join(out, 'base', 'default', cond, rep, f'{g}-genes', 'species.tre')
            if not os.path.isfile(seed):
                raise RuntimeError(f'no base default tree to seed from: {seed}')
            cmd += ['-c', seed]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        est = run_all.read_trees(os.path.join(d, 'species.tre'))[0]
        return {'variant': v, 'arm': arm, 'condition': cond, 'replicate': rep, 'genes': str(g),
                'fn': f'{fn_rate(run_all.read_trees(true_path)[0], est):.6f}',
                'top_changed': top_changed(os.path.join(d, 'subproblems.tsv'))}

    new = open(csv_path, 'a', newline='')
    w = csv.DictWriter(new, FIELDS)
    if os.path.getsize(csv_path) == 0:
        w.writeheader()
    done = 0
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for fut in cf.as_completed([ex.submit(run, j) for j in jobs]):
            try:
                r = fut.result()
            except (subprocess.CalledProcessError, RuntimeError) as e:
                print(f'FAILED: {e}', file=sys.stderr, flush=True)
                continue
            rows[(r['variant'], r['arm'], r['condition'], r['replicate'], r['genes'])] = r
            w.writerow(r)
            new.flush()
            done += 1
            if done % 200 == 0 or done == len(jobs):
                print(f'{done}/{len(jobs)} runs', flush=True)
    new.close()
    summary(list(rows.values()), os.path.join(out, 'summary.txt'))


if __name__ == '__main__':
    main()
