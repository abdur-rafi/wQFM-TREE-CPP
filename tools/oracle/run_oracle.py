"""Oracle experiment: would better initial bipartitions improve wQFM-TREE?

For every replicate, wqfm-tree runs twice on the same gene trees:
  consensus  the default: each subproblem starts from the best-scoring edge of the greedy consensus
             of the gene trees
  oracle     -c <true species tree>: each subproblem starts from the best-scoring edge of the TRUE
             species tree (projected onto the subproblem's taxa and dummy taxa); FM passes as usual
Both estimates and the true tree are then scored with
  quartet score  the fraction of gene-tree quartets a tree agrees with (build/qscore, exact; the
                 score ASTRAL reports)
  FN rate        the fraction of true-tree branches missing from an estimate
Oracle estimates with higher quartet scores or lower FN rates than the consensus ones mean that better
initial bipartitions would help. If the estimates already score at least as high as the true tree,
a better optimiser of the quartet score would not make them more accurate.

Each run also writes a subproblem log (wqfm-tree --subproblem-log): for every subproblem, the score
of its initial bipartition, the score after FM, and whether FM changed the bipartition. In the oracle
run the top-level subproblem starts from the best-scoring true-tree edge, so if FM changes it there,
the wQFM objective scores a bipartition that is not in the true tree higher than every true one.

usage: python3 -I tools/oracle/run_oracle.py -i INPUT -o OUTPUT [-t THREADS] [options]

  -i, --input DIR     Data in the run_all.py layout:
                      DIR/estimated-gene-trees/<condition>/<replicate>/estimatedgenetre
                      DIR/true_species_trees/<condition>/<replicate>/s_tree.trees
  -o, --output DIR    Results folder (created if needed), containing:
                      OUTPUT/results.csv    one row per replicate and gene count
                      OUTPUT/summary.txt    per condition and gene count
                      OUTPUT/setup.json     command, build and date of each invocation
                      OUTPUT/<condition>/<replicate>/<genes>-genes/{consensus,oracle}/species.tre, run.log,
                                                                         run.json, subproblems.tsv
                      OUTPUT/<condition>/<replicate>/<genes>-genes/qscore.tsv
  -t, --threads N     Threads per run (default: 1; the trees do not depend on it)
  --genes LIST        Gene counts; each run uses the first k gene trees (default: all)
  --conditions LIST   Conditions (comma-separated; default: the nine of run_all.py)
  --replicates LIST   Replicates, e.g. 1-10,15 (default: all numeric folders)
  --true-trees DIR    True species trees folder, if not DIR/true_species_trees
  --gene-trees-dir NAME  Gene-tree folder in DIR (default: estimated-gene-trees; e.g. true-gene-trees
                      with --gene-file truegenetrees for the simulated true gene trees)
  --gene-file NAME    Gene-tree filename (default: estimatedgenetre)
  --extra ARGS        Extra wqfm-tree arguments for both runs, as one string (--extra=--relative-eps)
  --heartbeat SECS    Echo a running run's latest log line every SECS (default: 60; 0: none)
  --overwrite         Rerun the selected jobs, replacing their outputs
  --summary-only      Only rebuild summary.txt from results.csv

Resuming: finished runs are kept and only missing steps run, unless --overwrite is given. Rows of
conditions and gene counts outside the selection stay in results.csv.

results.csv columns: condition, replicate, genes, taxa; qs_true, qs_consensus, qs_oracle (normalised
quartet scores); fn_consensus, fn_oracle, fp_consensus, fp_oracle; rf_between (normalised RF distance
between the two estimates); seconds_consensus, seconds_oracle; quartets (resolved gene-tree quartets)
and agree_true, agree_consensus, agree_oracle (exact counts behind the quartet scores);
fm_changed_<run> (fraction of subproblems whose bipartition FM changed), top_changed_<run> (1 if FM
changed the top-level bipartition) and top_gain_<run> (FM's relative score gain at the top level).
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import shlex
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HERE, '..', 'astral2'))
sys.path.insert(0, os.path.join(HERE, '..', 'compare'))
import rf  # noqa: E402
import run_all  # noqa: E402

sys.setrecursionlimit(100000)

ARMS = ('consensus', 'oracle')
FIELDS = ['condition', 'replicate', 'genes', 'taxa', 'qs_true', 'qs_consensus', 'qs_oracle',
          'fn_consensus', 'fn_oracle', 'fp_consensus', 'fp_oracle', 'rf_between',
          'seconds_consensus', 'seconds_oracle', 'quartets', 'agree_true', 'agree_consensus', 'agree_oracle',
          'fm_changed_consensus', 'fm_changed_oracle', 'top_changed_consensus', 'top_changed_oracle',
          'top_gain_consensus', 'top_gain_oracle']


def remove(*paths):
    for p in paths:
        try:
            os.remove(p)
        except FileNotFoundError:
            pass


def finished(run_dir):
    return all(os.path.isfile(os.path.join(run_dir, f)) for f in ('species.tre', 'run.json', 'subproblems.tsv'))


def fm_stats(path):
    """(fraction of subproblems whose bipartition FM changed, top level changed (0/1), FM's relative
    score gain at the top level) from a --subproblem-log file."""
    n = changed = top_changed = 0
    top_gain = 0.0
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            _, dummies, s0, s1, _, ch = line.split('\t')
            n += 1
            changed += int(ch)
            if int(dummies) == 0:                # the top-level subproblem
                top_changed = int(ch)
                top_gain = (float(s1) - float(s0)) / abs(float(s0)) if float(s0) else 0.0
    return (changed / n if n else 0.0), top_changed, top_gain


def qscore(binary, gene_trees, trees, threads):
    """(agree, resolved) for each species tree file."""
    out = subprocess.run([binary, '-t', str(threads), '-g', gene_trees] + trees, check=True,
                         capture_output=True, text=True).stdout
    rows = [line.split('\t') for line in out.strip().splitlines()]
    return [(int(r[0]), int(r[1])) for r in rows]


def splits_of(newick):
    return rf.splits(rf.parse(newick))


def row_for(cond, rep, genes, d, true_newick, scores):
    st, taxa = splits_of(true_newick)
    row = {'condition': cond, 'replicate': rep, 'genes': genes, 'taxa': len(taxa)}
    est = {}
    for arm in ARMS:
        newick = run_all.read_trees(os.path.join(d, arm, 'species.tre'))[0]
        se, te = splits_of(newick)
        if te != taxa:
            raise RuntimeError(f'{d}/{arm}: estimate and true tree have different taxa')
        est[arm] = se
        row[f'fn_{arm}'] = f'{len(st - se) / len(st):.6f}'
        row[f'fp_{arm}'] = f'{len(se - st) / len(se):.6f}'
        row[f'seconds_{arm}'] = json.load(open(os.path.join(d, arm, 'run.json')))['seconds']
        ch, top_ch, top_gain = fm_stats(os.path.join(d, arm, 'subproblems.tsv'))
        row[f'fm_changed_{arm}'] = f'{ch:.6f}'
        row[f'top_changed_{arm}'] = top_ch
        row[f'top_gain_{arm}'] = f'{top_gain:.3e}'
    row['rf_between'] = f'{len(est["consensus"] ^ est["oracle"]) / (2 * (len(taxa) - 3)):.6f}'
    resolved = scores['true'][1]
    row['quartets'] = resolved
    for name in ('true',) + ARMS:
        agree, res = scores[name]
        if res != resolved:
            raise RuntimeError(f'{d}: quartet totals differ between trees')
        row[f'agree_{name}'] = agree
        row[f'qs_{name}'] = f'{agree / resolved:.8f}'
    return row


def read_qscore(path):
    scores = {}
    with open(path) as f:
        for line in f:
            name, agree, res = line.split('\t')
            scores[name] = (int(agree), int(res))
    return scores


def sign_test(wins, losses):
    """Two-sided exact sign test p-value (ties dropped)."""
    n = wins + losses
    if n == 0:
        return 1.0
    k = min(wins, losses)
    p = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * p)


def summary(rows, path):
    groups = {}
    for r in rows:
        groups.setdefault((r['condition'], int(r['genes'])), []).append(r)
    head = (f'{"condition":30s} {"genes":>5s} {"runs":>4s} | {"QS true":>8s} {"QS cons":>8s} {"QS orcl":>8s} '
            f'{"orcl>":>5s} {"=":>3s} {"<":>3s} {"p":>7s} {"cons>=true":>10s} | {"FN cons":>7s} {"FN orcl":>7s} '
            f'{"better":>6s} {"=":>3s} {"worse":>5s} {"p":>7s} | {"RF c-o":>6s} {"same":>4s}')
    lines = [head, '-' * len(head)]
    for (cond, genes), rs in sorted(groups.items(), key=lambda x: (int(x[0][0].split('.')[1]), x[0][0], x[0][1])):
        mean = lambda col: statistics.mean(float(r[col]) for r in rs)
        qs_w = sum(int(r['agree_oracle']) > int(r['agree_consensus']) for r in rs)
        qs_l = sum(int(r['agree_oracle']) < int(r['agree_consensus']) for r in rs)
        fn_w = sum(float(r['fn_oracle']) < float(r['fn_consensus']) for r in rs)
        fn_l = sum(float(r['fn_oracle']) > float(r['fn_consensus']) for r in rs)
        beats = sum(int(r['agree_consensus']) >= int(r['agree_true']) for r in rs)
        same = sum(float(r['rf_between']) == 0 for r in rs)
        lines.append(f'{cond:30s} {genes:5d} {len(rs):4d} | {mean("qs_true"):8.5f} {mean("qs_consensus"):8.5f} '
                     f'{mean("qs_oracle"):8.5f} {qs_w:5d} {len(rs) - qs_w - qs_l:3d} {qs_l:3d} '
                     f'{sign_test(qs_w, qs_l):7.1e} {beats:10d} | {mean("fn_consensus"):7.4f} {mean("fn_oracle"):7.4f} '
                     f'{fn_w:6d} {len(rs) - fn_w - fn_l:3d} {fn_l:5d} {sign_test(fn_w, fn_l):7.1e} | '
                     f'{mean("rf_between"):6.4f} {same:4d}')
    head2 = (f'{"condition":30s} {"genes":>5s} {"runs":>4s} | {"FM changed (% subproblems)":>26s} | '
             f'{"top level changed by FM":>23s} | {"top-level FM gain":>17s}')
    lines += ['', head2, f'{"":30s} {"":5s} {"":4s} | {"cons":>12s} {"orcl":>13s} | {"cons":>11s} {"orcl":>11s} | '
              f'{"orcl, mean":>17s}', '-' * len(head2)]
    for (cond, genes), rs in sorted(groups.items(), key=lambda x: (int(x[0][0].split('.')[1]), x[0][0], x[0][1])):
        mean = lambda col: statistics.mean(float(r[col]) for r in rs)
        lines.append(f'{cond:30s} {genes:5d} {len(rs):4d} | {100 * mean("fm_changed_consensus"):11.1f}% '
                     f'{100 * mean("fm_changed_oracle"):12.1f}% | {sum(int(r["top_changed_consensus"]) for r in rs):11d} '
                     f'{sum(int(r["top_changed_oracle"]) for r in rs):11d} | {mean("top_gain_oracle"):17.2e}')
    lines += ['',
              'QS: normalised quartet score (fraction of gene-tree quartets the tree agrees with); cons = default',
              '(greedy-consensus initial bipartitions), orcl = initial bipartitions from the true species tree.',
              'orcl> / = / <: replicates where the oracle estimate has a higher / equal / lower quartet score;',
              'p: two-sided exact sign test (ties dropped). cons>=true: replicates where the default estimate',
              'scores at least as high as the true tree. FN better / = / worse: oracle FN rate lower / equal /',
              'higher. RF c-o: mean normalised RF distance between the two estimates; same: identical trees.',
              'FM changed: subproblems where FM moved away from the initial bipartition. Top level changed: runs',
              'where it did so in the top-level subproblem; for the oracle runs, the objective then scores a',
              'bipartition that is not in the true tree above every true-tree edge. Top-level FM gain: FM\'s',
              'relative score gain there.']
    text = '\n'.join(lines) + '\n'
    run_all.write_atomic(path, text)
    print(text)


def write_csv(rows, path):
    ordered = sorted(rows.values(), key=lambda r: (r['condition'], r['replicate'], int(r['genes'])))
    tmp = path + '.part'
    with open(tmp, 'w', newline='') as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        w.writerows(ordered)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def setup_record(out, args):
    """Append how this invocation was run (command, repository state, date) to OUTPUT/setup.json."""
    def git(*cmd):
        try:
            return subprocess.run(['git', '-C', REPO] + list(cmd), capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return ''
    rec = {'date': datetime.datetime.now().isoformat(timespec='seconds'),
           'command': ' '.join(shlex.quote(x) for x in [sys.executable] + sys.argv),
           'commit': git('rev-parse', 'HEAD'),
           'uncommitted_changes': git('status', '--porcelain', '--untracked-files=no') != '',
           'binary': args.binary, 'qscore': args.qscore, 'threads': args.threads, 'extra': args.extra}
    for key in ('binary', 'qscore'):             # pins the exact build, committed or not
        with open(getattr(args, key), 'rb') as f:
            rec[f'{key}_sha256'] = hashlib.sha256(f.read()).hexdigest()
    path = os.path.join(out, 'setup.json')
    hist = json.load(open(path)) if os.path.isfile(path) else []
    hist.append(rec)
    run_all.write_atomic(path, json.dumps(hist, indent=1) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-i', '--input')
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('-t', '--threads', type=int, default=1)
    ap.add_argument('--true-trees')
    ap.add_argument('--genes')
    ap.add_argument('--conditions')
    ap.add_argument('--replicates')
    ap.add_argument('--gene-trees-dir', default=run_all.GENE_TREES_DIR)
    ap.add_argument('--gene-file', default='estimatedgenetre')
    ap.add_argument('--binary', default=os.path.join(REPO, 'build', 'wqfm-tree'))
    ap.add_argument('--qscore', default=os.path.join(REPO, 'build', 'qscore'))
    ap.add_argument('--extra', default='')
    ap.add_argument('--heartbeat', type=float, default=60)
    ap.add_argument('--overwrite', action='store_true')
    ap.add_argument('--summary-only', action='store_true')
    a = ap.parse_args()
    out = os.path.abspath(a.output)
    csv_path = os.path.join(out, 'results.csv')
    if a.summary_only:
        with open(csv_path, newline='') as f:
            summary(list(csv.DictReader(f)), os.path.join(out, 'summary.txt'))
        return
    if not a.input:
        ap.error('-i/--input is required')
    data = os.path.abspath(a.input)
    inp = os.path.join(data, a.gene_trees_dir)
    true_root = os.path.abspath(a.true_trees) if a.true_trees else os.path.join(data, run_all.TRUE_TREES_DIR)
    for exe in (a.binary, a.qscore):
        if not os.path.isfile(exe):
            sys.exit(f'no executable at {exe} (run make)')
    if not os.path.isdir(inp):
        sys.exit(f'no {a.gene_trees_dir} folder in {data}')
    if not os.path.isdir(true_root):
        sys.exit(f'no true species trees at {true_root}')
    os.makedirs(out, exist_ok=True)
    setup_record(out, a)

    conds = a.conditions.split(',') if a.conditions else run_all.CONDITIONS
    for cond in conds:
        if cond not in run_all.CONDITIONS:
            sys.exit(f'unknown condition {cond}; expected one of:\n  ' + '\n  '.join(run_all.CONDITIONS))
    want_reps = {f'{r:02d}' for r in run_all.parse_list(a.replicates)} if a.replicates else None
    jobs = []
    for cond in conds:
        if not os.path.isdir(os.path.join(inp, cond)):
            print(f'skip {cond}: no folder {os.path.join(inp, cond)}', file=sys.stderr)
            continue
        for rep in run_all.replicates(os.path.join(inp, cond)):
            if want_reps is not None and rep not in want_reps:
                continue
            src = os.path.join(inp, cond, rep, a.gene_file)
            true_path = os.path.join(true_root, cond, rep, 's_tree.trees')
            if not os.path.isfile(src) or not os.path.isfile(true_path):
                print(f'skip {cond}/{rep}: missing {src if not os.path.isfile(src) else true_path}', file=sys.stderr)
                continue
            n = len(run_all.read_trees(src))
            for genes in (run_all.parse_list(a.genes) if a.genes else [n]):
                if genes > n:
                    print(f'skip {cond}/{rep}: {genes} genes requested, file has {n}', file=sys.stderr)
                    continue
                jobs.append((cond, rep, genes, n, src, true_path))

    def job_dir(cond, rep, genes):
        return os.path.join(out, cond, rep, f'{genes}-genes')

    rows = {}                                    # rows of other conditions or gene counts are kept
    if os.path.exists(csv_path):
        with open(csv_path, newline='') as f:
            for r in csv.DictReader(f):
                if all(r.get(k) not in (None, '') for k in FIELDS):
                    rows[(r['condition'], r['replicate'], int(r['genes']))] = r
    todo = []
    for job in jobs:
        cond, rep, genes, n, src, true_path = job
        d = job_dir(cond, rep, genes)
        if a.overwrite:                          # invalidate up front: an interrupted sweep reruns them
            for arm in ARMS:
                remove(os.path.join(d, arm, 'species.tre'), os.path.join(d, arm, 'run.json'))
            remove(os.path.join(d, 'qscore.tsv'))
        if all(finished(os.path.join(d, arm)) for arm in ARMS) and os.path.isfile(os.path.join(d, 'qscore.tsv')):
            true_newick = run_all.read_trees(true_path)[0]
            rows[(cond, rep, genes)] = row_for(cond, rep, genes, d, true_newick, read_qscore(os.path.join(d, 'qscore.tsv')))
        else:
            rows.pop((cond, rep, genes), None)
            todo.append(job)
    write_csv(rows, csv_path)
    print(f'{len(jobs)} jobs: {len(jobs) - len(todo)} already finished, {len(todo)} to go; '
          f'{a.threads} thread(s) per run' + (' [overwrite]' if a.overwrite else ''), flush=True)

    extra = a.extra.split()
    t_start = time.time()
    for k, (cond, rep, genes, n, src, true_path) in enumerate(todo, 1):
        d = job_dir(cond, rep, genes)
        os.makedirs(d, exist_ok=True)
        gt = src
        if genes < n:                            # the first `genes` gene trees
            gt = os.path.join(d, 'gene-trees.tre')
            if not os.path.isfile(gt):
                run_all.write_atomic(gt, '\n'.join(run_all.read_trees(src)[:genes]) + '\n')
        label = f'[{len(jobs) - len(todo) + k}/{len(jobs)}] {cond}/{rep}, {genes} genes'
        print(f'{label}: starting ({time.strftime("%H:%M:%S")})', flush=True)
        try:
            for arm in ARMS:
                rd = os.path.join(d, arm)
                if finished(rd):
                    continue
                remove(os.path.join(d, 'qscore.tsv'))
                args = extra + ['--subproblem-log', os.path.join(rd, 'subproblems.tsv')]
                run_all.run_one(a.binary, gt, rd, a.threads, args + (['-c', true_path] if arm == 'oracle' else []),
                                a.heartbeat)
            qs_path = os.path.join(d, 'qscore.tsv')
            if not os.path.isfile(qs_path):
                names = ('true',) + ARMS
                res = qscore(a.qscore, gt, [true_path] + [os.path.join(d, arm, 'species.tre') for arm in ARMS], a.threads)
                run_all.write_atomic(qs_path, ''.join(f'{nm}\t{ag}\t{rs}\n' for nm, (ag, rs) in zip(names, res)))
            row = row_for(cond, rep, genes, d, run_all.read_trees(true_path)[0], read_qscore(qs_path))
        except (RuntimeError, subprocess.CalledProcessError) as e:
            print(f'{label}: FAILED: {e}', file=sys.stderr, flush=True)
            continue
        rows[(cond, rep, genes)] = row
        write_csv(rows, csv_path)
        elapsed = time.time() - t_start
        print(f'{label}: QS true {float(row["qs_true"]):.5f}, consensus {float(row["qs_consensus"]):.5f}, '
              f'oracle {float(row["qs_oracle"]):.5f}; FN consensus {row["fn_consensus"]}, oracle {row["fn_oracle"]} '
              f'| {elapsed / 3600:.2f} h elapsed, ~{elapsed / k * (len(todo) - k) / 3600:.2f} h left', flush=True)
    if rows:
        summary(rows.values(), os.path.join(out, 'summary.txt'))


if __name__ == '__main__':
    main()
