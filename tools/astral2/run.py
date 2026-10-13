"""Run the C++ wQFM-TREE on the ASTRAL-II simulated datasets and score the results.

usage: python3 -I tools/astral2/run.py --data DIR --out DIR [options]

  --data DIR        extracted estimated-gene-trees archive (DIR/<condition>/<rep>/estimatedgenetre)
  --out DIR         results: per-run species trees and logs, results.csv, summary
  --results DIR     the repo's results/Astral_2_datasets (true species trees, published trees);
                    default: ../results/Astral_2_datasets relative to this repo
  --conds LIST      comma-separated conditions (default: the eight analysed in the paper)
  --reps LIST       e.g. 1-5,7 (default: the replicates with published trees)
  --genes LIST      gene counts, the first k trees of each file (default 50,200,1000)
  --threads N       passed to -t (default: binary's default)
  --input NAME      gene-tree file name (default estimatedgenetre; e.g. estimatedgenetre.halfresolved)
  --extra ARGS      extra arguments for the binary, as one string
  --force           rerun runs whose species tree exists

Accuracy is the missing-branch rate against the true species tree (the FN rate of the published
score files). The published wQFM-TREE trees used PAUP's greedy consensus; this uses the built-in one.
"""
import argparse
import csv
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HERE, '..', 'compare'))
import rf  # noqa: E402

sys.setrecursionlimit(100000)

PAPER_CONDS = ['model.200.500000.0.000001', 'model.200.500000.0.0000001',
               'model.200.2000000.0.000001', 'model.200.2000000.0.0000001',
               'model.200.10000000.0.000001', 'model.200.10000000.0.0000001',
               'model.500.2000000.0.000001', 'model.1000.2000000.0.000001']
FIELDS = ['cond', 'rep', 'genes', 'threads', 'seconds', 'max_rss_mb', 'fn_true', 'pub_fn_true', 'rf_pub']


def published_name(genes):
    return 'ewqfm-paup-nq-A-tree.txt' if genes == 1000 else f'ewqfm-{genes}-tree.txt'


def parse_list(s):
    out = []
    for part in s.split(','):
        if '-' in part:
            a, b = map(int, part.split('-'))
            out.extend(range(a, b + 1))
        else:
            out.append(int(part))
    return out


def splits_of(path):
    if not os.path.exists(path):
        return None
    return rf.splits(rf.parse(rf.first_tree(path)))


def fn_rate(true, est):
    (st, tt), (se, te) = true, est
    if tt != te:
        raise ValueError('different taxon sets')
    return len(st - se) / len(st)


def rf_norm(a, b):
    (sa, ta), (sb, tb) = a, b
    if ta != tb:
        raise ValueError('different taxon sets')
    return len(sa ^ sb) / (2 * (len(ta) - 3))


def gene_file(src, genes, dst):
    """Path of a file holding the first `genes` trees of src (src itself if it has exactly that many)."""
    with open(src) as f:
        lines = [line for line in f if '(' in line]
    if len(lines) < genes:
        raise ValueError(f'{src}: only {len(lines)} trees')
    if len(lines) == genes:
        return src
    if not os.path.exists(dst):
        with open(dst + '.part', 'w') as f:
            f.writelines(lines[:genes])
        os.replace(dst + '.part', dst)
    return dst


def run(binary, gt, out_tree, log, threads, extra):
    cmd = [binary, '-i', gt, '-o', out_tree] + (['-t', str(threads)] if threads else []) + extra
    with open(log, 'w') as lf:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
        secs = time.perf_counter() - t0
    p.returncode = os.waitstatus_to_exitcode(status)
    if p.returncode != 0:
        raise RuntimeError(f'{" ".join(cmd)} failed ({p.returncode}), see {log}')
    rss = ru.ru_maxrss / (1 << 20) if sys.platform == 'darwin' else ru.ru_maxrss / 1024
    return secs, rss


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--results', default=os.path.join(REPO, '..', 'results', 'Astral_2_datasets'))
    ap.add_argument('--conds', default=','.join(PAPER_CONDS))
    ap.add_argument('--reps')
    ap.add_argument('--genes', default='50,200,1000')
    ap.add_argument('--threads', type=int, default=0)
    ap.add_argument('--input', default='estimatedgenetre')
    ap.add_argument('--extra', default='')
    ap.add_argument('--binary', default=os.path.join(REPO, 'build', 'wqfm-tree'))
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    csv_path = os.path.join(a.out, 'results.csv')
    rows = {}
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            for r in csv.DictReader(f):
                rows[(r['cond'], r['rep'], int(r['genes']))] = r
    extra = a.extra.split()

    for cond in a.conds.split(','):
        pub_dir = os.path.join(a.results, 'outputs', cond)
        if a.reps:
            reps = [f'{r:02d}' for r in parse_list(a.reps)]
        else:
            reps = sorted(d for d in os.listdir(pub_dir) if d.isdigit())
        for rep in reps:
            src = os.path.join(a.data, cond, rep, a.input)
            if not os.path.exists(src):
                print(f'skip {cond}/{rep}: no {src}', file=sys.stderr)
                continue
            true = splits_of(os.path.join(a.results, 'true-species-trees', cond, rep, 'sp-cleaned'))
            od = os.path.join(a.out, cond, rep)
            os.makedirs(od, exist_ok=True)
            for genes in parse_list(a.genes):
                key = (cond, rep, genes)
                out_tree = os.path.join(od, f'species-{genes}.tre')
                if key in rows and os.path.exists(out_tree) and not a.force:
                    continue
                gt = gene_file(src, genes, os.path.join(od, f'gt-{genes}.tre'))
                secs, rss = run(a.binary, gt, out_tree, os.path.join(od, f'run-{genes}.log'), a.threads, extra)
                est = splits_of(out_tree)
                pub = splits_of(os.path.join(pub_dir, rep, published_name(genes)))
                row = {'cond': cond, 'rep': rep, 'genes': genes, 'threads': a.threads or '',
                       'seconds': f'{secs:.3f}', 'max_rss_mb': f'{rss:.1f}',
                       'fn_true': f'{fn_rate(true, est):.6f}' if true else '',
                       'pub_fn_true': f'{fn_rate(true, pub):.6f}' if true and pub else '',
                       'rf_pub': f'{rf_norm(est, pub):.6f}' if pub else ''}
                rows[key] = row
                print(' '.join(f'{k}={row[k]}' for k in FIELDS), flush=True)
                with open(csv_path + '.part', 'w', newline='') as f:
                    w = csv.DictWriter(f, FIELDS)
                    w.writeheader()
                    w.writerows(sorted(rows.values(), key=lambda r: (r['cond'], r['rep'], int(r['genes']))))
                os.replace(csv_path + '.part', csv_path)

    summary(rows.values(), os.path.join(a.out, 'summary.txt'))


def summary(rows, path):
    groups = {}
    for r in rows:
        groups.setdefault((r['cond'], int(r['genes'])), []).append(r)
    lines = ['condition                      genes reps  FN C++  FN published  same as published  '
             'mean s   max s  max RSS MB']
    for (cond, genes), rs in sorted(groups.items()):
        both = [r for r in rs if r['fn_true'] and r['pub_fn_true']]
        fn = statistics.mean(float(r['fn_true']) for r in both) if both else float('nan')
        pfn = statistics.mean(float(r['pub_fn_true']) for r in both) if both else float('nan')
        same = sum(1 for r in rs if r['rf_pub'] and float(r['rf_pub']) == 0)
        npub = sum(1 for r in rs if r['rf_pub'])
        secs = [float(r['seconds']) for r in rs]
        lines.append(f'{cond:30s} {genes:5d} {len(rs):4d}  {fn:6.4f}  {pfn:12.4f}  {same:8d}/{npub:<8d}  '
                     f'{statistics.mean(secs):6.2f}  {max(secs):6.2f}  {max(float(r["max_rss_mb"]) for r in rs):10.0f}')
    text = '\n'.join(lines) + '\n'
    with open(path, 'w') as f:
        f.write(text)
    print(text)


if __name__ == '__main__':
    main()
