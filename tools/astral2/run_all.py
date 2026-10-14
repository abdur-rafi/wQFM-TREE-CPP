"""Run wQFM-TREE (C++) on every model condition and replicate of the ASTRAL-II simulated data, resumably.

usage: python3 -I tools/astral2/run_all.py -i INPUT -o OUTPUT [-t THREADS] [options]

  -i, --input DIR     estimated gene trees: DIR/<condition>/<replicate>/estimatedgenetre
  -o, --output DIR    results (created if needed), laid out as
                        OUTPUT/results.csv                       one row per finished run
                        OUTPUT/summary.txt                       per condition and gene count
                        OUTPUT/<condition>/<replicate>/<genes>-genes/species.tre   estimated species tree
                                                                 /run.log        the program's output
                                                                 /run.json       time and memory
  -t, --threads N     threads per run (default 1)
  --true-trees DIR    true species trees: DIR/<condition>/<replicate>/s_tree.trees
                      (default: the folder true-specis-trees next to INPUT)
  --genes LIST        gene counts, each run on the first k gene trees of the file (default: all of them)
  --conditions LIST   only these conditions (comma-separated; default: every folder in INPUT)
  --replicates LIST   only these replicates, e.g. 1-10,15 (default: every folder)
  --gene-file NAME    gene-tree file name (default estimatedgenetre)
  --binary PATH       the wqfm-tree executable (default: build/wqfm-tree of this repository)
  --extra ARGS        more arguments for wqfm-tree, as one string (e.g. "--relative-eps")

Resuming: run the same command again. Finished runs (species.tre and run.json present) are skipped,
an interrupted run is started over, and results.csv is repaired from the finished runs, so a
shutdown at any moment loses at most the run in progress.

results.csv columns: condition, replicate, genes, taxa, threads, seconds, max_rss_mb,
fn_rate (true-tree branches missing from the estimate), fp_rate (estimated branches not in the
true tree), rf (normalised Robinson-Foulds distance), estimated_tree, true_tree (Newick).
"""
import argparse
import csv
import json
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
csv.field_size_limit(sys.maxsize)

FIELDS = ['condition', 'replicate', 'genes', 'taxa', 'threads', 'seconds', 'max_rss_mb',
          'fn_rate', 'fp_rate', 'rf', 'estimated_tree', 'true_tree']


def parse_list(s):
    out = []
    for part in s.split(','):
        if '-' in part:
            a, b = map(int, part.split('-'))
            out.extend(range(a, b + 1))
        elif part:
            out.append(int(part))
    return out


def subdirs(path):
    return sorted(d for d in os.listdir(path) if not d.startswith('.') and os.path.isdir(os.path.join(path, d)))


def read_trees(path):
    with open(path) as f:
        return [line.strip() for line in f if '(' in line]


def write_atomic(path, text):
    tmp = path + '.part'
    with open(tmp, 'w') as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def compare(true_newick, est_newick):
    """(taxa, fn_rate, fp_rate, rf) of the estimate against the true tree; rates blank if the taxa differ."""
    st, tt = rf.splits(rf.parse(true_newick))
    se, te = rf.splits(rf.parse(est_newick))
    if tt != te:
        return len(te), '', '', ''
    fn = len(st - se) / len(st) if st else 0.0
    fp = len(se - st) / len(se) if se else 0.0
    return len(te), f'{fn:.6f}', f'{fp:.6f}', f'{len(st ^ se) / (2 * (len(tt) - 3)):.6f}'


def run_one(binary, gene_trees, run_dir, threads, extra):
    """Run wqfm-tree in run_dir; species.tre appears only when the run has finished."""
    os.makedirs(run_dir, exist_ok=True)
    tree, tmp = os.path.join(run_dir, 'species.tre'), os.path.join(run_dir, 'species.tre.part')
    cmd = [binary, '-i', gene_trees, '-o', tmp, '-t', str(threads)] + extra
    with open(os.path.join(run_dir, 'run.log'), 'w') as log:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(p.pid, 0)
        secs = time.perf_counter() - t0
    code = os.waitstatus_to_exitcode(status)
    if code != 0:
        raise RuntimeError(f'{" ".join(cmd)} exited with {code}, see {run_dir}/run.log')
    rss = usage.ru_maxrss / (1 << 20) if sys.platform == 'darwin' else usage.ru_maxrss / 1024
    write_atomic(os.path.join(run_dir, 'run.json'), json.dumps({'seconds': round(secs, 3), 'max_rss_mb': round(rss, 1), 'threads': threads}))
    os.replace(tmp, tree)


def row_for(cond, rep, genes, run_dir, true_newick):
    est = read_trees(os.path.join(run_dir, 'species.tre'))[0]
    meta = json.load(open(os.path.join(run_dir, 'run.json')))
    taxa, fn, fp, d = compare(true_newick, est) if true_newick else (len(rf.splits(rf.parse(est))[1]), '', '', '')
    return {'condition': cond, 'replicate': rep, 'genes': genes, 'taxa': taxa, 'threads': meta['threads'],
            'seconds': meta['seconds'], 'max_rss_mb': meta['max_rss_mb'], 'fn_rate': fn, 'fp_rate': fp, 'rf': d,
            'estimated_tree': est, 'true_tree': true_newick or ''}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('-t', '--threads', type=int, default=1)
    ap.add_argument('--true-trees')
    ap.add_argument('--genes')
    ap.add_argument('--conditions')
    ap.add_argument('--replicates')
    ap.add_argument('--gene-file', default='estimatedgenetre')
    ap.add_argument('--binary', default=os.path.join(REPO, 'build', 'wqfm-tree'))
    ap.add_argument('--extra', default='')
    a = ap.parse_args()

    inp, out = os.path.abspath(a.input), os.path.abspath(a.output)
    true_root = a.true_trees or os.path.join(os.path.dirname(inp), 'true-specis-trees')
    if not os.path.isfile(a.binary):
        sys.exit(f'no wqfm-tree executable at {a.binary} (run make, or pass --binary)')
    if not os.path.isdir(true_root):
        print(f'warning: no true trees at {true_root}; accuracy columns stay empty', file=sys.stderr)
    os.makedirs(out, exist_ok=True)

    # every (condition, replicate, genes) to run, in order
    conds = a.conditions.split(',') if a.conditions else subdirs(inp)
    want_reps = {f'{r:02d}' for r in parse_list(a.replicates)} if a.replicates else None
    jobs = []
    for cond in conds:
        for rep in subdirs(os.path.join(inp, cond)):
            if want_reps is not None and rep not in want_reps:
                continue
            src = os.path.join(inp, cond, rep, a.gene_file)
            if not os.path.isfile(src):
                print(f'skip {cond}/{rep}: no {a.gene_file}', file=sys.stderr)
                continue
            for genes in (parse_list(a.genes) if a.genes else [None]):
                jobs.append((cond, rep, genes, src))

    # finished runs, and a results.csv rebuilt from them (repairs a row cut off by a shutdown)
    csv_path = os.path.join(out, 'results.csv')
    def complete(r):                     # a row cut off by a shutdown fails one of these
        if None in r or any(r.get(k) is None for k in FIELDS) or None in r.values():
            return False
        try:
            float(r['seconds']), float(r['max_rss_mb']), int(r['genes'])
        except ValueError:
            return False
        return r['estimated_tree'].endswith(';') and (r['true_tree'] == '' or r['true_tree'].endswith(';'))

    rows = {}
    if os.path.exists(csv_path):
        with open(csv_path, newline='') as f:
            for r in csv.DictReader(f):
                if complete(r):
                    rows[(r['condition'], r['replicate'], r['genes'])] = r
    true_cache = {}

    def true_tree(cond, rep):
        if (cond, rep) not in true_cache:
            p = os.path.join(true_root, cond, rep, 's_tree.trees')
            true_cache[(cond, rep)] = read_trees(p)[0] if os.path.isfile(p) else None
        return true_cache[(cond, rep)]

    def run_dir(cond, rep, genes, n):
        return os.path.join(out, cond, rep, f'{genes or n}-genes')

    gene_counts = {}

    def gene_count(src):
        if src not in gene_counts:
            gene_counts[src] = len(read_trees(src))
        return gene_counts[src]

    todo = []
    for cond, rep, genes, src in jobs:
        n = gene_count(src)
        if genes is not None and genes > n:
            print(f'skip {cond}/{rep}: {genes} genes requested, file has {n}', file=sys.stderr)
            continue
        key = (cond, rep, str(genes or n))
        d = run_dir(cond, rep, genes, n)
        finished = os.path.isfile(os.path.join(d, 'species.tre')) and os.path.isfile(os.path.join(d, 'run.json'))
        if finished and key not in rows:
            rows[key] = row_for(cond, rep, genes or n, d, true_tree(cond, rep))
        if not finished:
            rows.pop(key, None)
            todo.append((cond, rep, genes, src, n))

    def save_csv():
        ordered = sorted(rows.values(), key=lambda r: (r['condition'], r['replicate'], int(r['genes'])))
        tmp = csv_path + '.part'
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, FIELDS)
            w.writeheader()
            w.writerows(ordered)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, csv_path)

    save_csv()
    done_before = len(jobs) - len(todo)
    print(f'{len(jobs)} runs: {done_before} already finished, {len(todo)} to go; {a.threads} thread(s) each', flush=True)

    # append new rows as runs finish (fsync'd), so the CSV never needs a full rewrite during the sweep
    extra = a.extra.split()
    t_start = time.time()
    with open(csv_path, 'a', newline='') as f:
        w = csv.DictWriter(f, FIELDS)
        for k, (cond, rep, genes, src, n) in enumerate(todo, 1):
            d = run_dir(cond, rep, genes, n)
            gt = src
            if genes is not None and genes < n:              # the first `genes` trees of the file
                gt = os.path.join(d, 'gene-trees.tre')
                if not os.path.isfile(gt):
                    os.makedirs(d, exist_ok=True)
                    write_atomic(gt, '\n'.join(read_trees(src)[:genes]) + '\n')
            try:
                run_one(a.binary, gt, d, a.threads, extra)
            except RuntimeError as e:
                print(f'FAILED {cond}/{rep} ({genes or n} genes): {e}', file=sys.stderr, flush=True)
                continue
            row = row_for(cond, rep, genes or n, d, true_tree(cond, rep))
            rows[(cond, rep, str(genes or n))] = row
            w.writerow(row)
            f.flush()
            os.fsync(f.fileno())
            elapsed = time.time() - t_start
            eta = elapsed / k * (len(todo) - k)
            print(f'[{done_before + k}/{len(jobs)}] {cond}/{rep} {genes or n} genes: {row["seconds"]} s, '
                  f'FN {row["fn_rate"] or "-"}, RF {row["rf"] or "-"} | ETA {eta / 3600:.1f} h', flush=True)

    save_csv()                                               # sorted
    summary(rows.values(), os.path.join(out, 'summary.txt'))


def summary(rows, path):
    groups = {}
    for r in rows:
        groups.setdefault((r['condition'], int(r['genes'])), []).append(r)
    lines = [f'{"condition":30s} {"genes":>5s} {"runs":>4s} {"mean FN":>8s} {"mean FP":>8s} {"mean RF":>8s} '
             f'{"mean s":>8s} {"max s":>8s} {"max MB":>7s}']
    for (cond, genes), rs in sorted(groups.items()):
        def mean(col):
            v = [float(r[col]) for r in rs if r[col] != '']
            return f'{statistics.mean(v):8.4f}' if v else f'{"-":>8s}'
        secs = [float(r['seconds']) for r in rs]
        lines.append(f'{cond:30s} {genes:5d} {len(rs):4d} {mean("fn_rate")} {mean("fp_rate")} {mean("rf")} '
                     f'{statistics.mean(secs):8.2f} {max(secs):8.2f} {max(float(r["max_rss_mb"]) for r in rs):7.0f}')
    text = '\n'.join(lines) + '\n'
    write_atomic(path, text)
    print(text)


if __name__ == '__main__':
    main()
