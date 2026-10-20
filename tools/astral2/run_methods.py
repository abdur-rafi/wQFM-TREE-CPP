"""Run weighted ASTRAL, TREE-QMC or wQFM-TREE on the ASTRAL-II simulated data, resumably, several runs at a time.

usage: python3 -I tools/astral2/run_methods.py -i INPUT -o OUTPUT --method METHOD --binary PATH [options]

  -i, --input DIR       Data folder, laid out as for run_all.py:
                        DIR/estimated-gene-trees/<condition>/<replicate>/estimatedgenetre
                        DIR/true_species_trees/<condition>/<replicate>/s_tree.trees
  -o, --output DIR      Results folder for ONE method and setting (created if needed), containing:
                        OUTPUT/results.csv                     one row per finished run
                        OUTPUT/summary.txt                     per condition and gene count
                        OUTPUT/settings.json                   method and settings of the folder
                        OUTPUT/<condition>/<replicate>/<genes>-genes/species.tre
                                                            /run.log
                                                            /run.json
  --method NAME         wastral    weighted ASTRAL (ASTER's bin/wastral)
                        tree-qmc   TREE-QMC (build/tree-qmc)
                        wqfm-tree  wQFM-TREE (C++, build/wqfm-tree), for comparisons under the same conditions
  --binary PATH         The method's executable
  --weight MODE         hybrid (default), support, length, or none
  -j, --jobs N          Runs at a time (default: 1). Every run uses a single thread.
  --genes LIST          Gene counts; each run uses the first k gene trees (default: all)
  --conditions LIST     Conditions to run (comma-separated; default: the nine of run_all.py)
  --replicates LIST     Replicates, e.g. 1-10,15 (default: all numeric folders)
  --gene-file NAME      Gene-tree filename (default: estimatedgenetre)
  --true-trees DIR      True species trees folder, if not DIR/true_species_trees
  --support-scale MIN MAX
                        Range of the gene-tree branch supports, mapped to 0..1 (default: 0 1, the
                        SH-like supports of the ASTRAL-II gene trees). Same for every method.
  --support-missing S   Support of a branch that has none, on the same scale (default: MAX, i.e.
                        fully supported, as in wQFM-TREE). In the ASTRAL-II gene trees these are the
                        branches above cherries of identical sequences, about 0.4 per tree.
  --timeout SECS        Stop a run after SECS seconds and count it as failed (default: no limit)
  --extra ARGS          Extra arguments for the method, as one string (recorded in settings.json)
  --heartbeat SECS      Echo the latest log line of each running run every SECS (default: 60; 0: none)
  --overwrite           Rerun selected jobs, replacing existing outputs and CSV records

Example: the comparison on the six 200-taxon conditions, 50 replicates, 1000 genes, 8 runs at a time
(bash or zsh, from the repository root; DATA holds estimated-gene-trees/ and true_species_trees/):
  DATA=/path/to/data
  C=model.200.500000.0.0000001,model.200.500000.0.000001,model.200.2000000.0.0000001,model.200.2000000.0.000001,model.200.10000000.0.0000001,model.200.10000000.0.000001
  S=(tools/astral2/run_methods.py -i "$DATA" --conditions "$C" --replicates 1-50 --genes 1000 -j 8)
  python3 -I "${S[@]}" -o runs/wastral-hybrid --method wastral   --binary /path/to/ASTER/bin/wastral
  python3 -I "${S[@]}" -o runs/treeqmc-hybrid --method tree-qmc  --binary /path/to/TREE-QMC/build/tree-qmc
  python3 -I "${S[@]}" -o runs/wqfm-hybrid    --method wqfm-tree --binary build/wqfm-tree
  python3 -I "${S[@]}" -o runs/wqfm-plain     --method wqfm-tree --binary build/wqfm-tree --weight none
  python3 -I tools/weighted/compare_fn.py runs/wastral-hybrid/results.csv runs/wqfm-hybrid/results.csv --labels wastral-h,wqfm-h
Rerun the same command to resume after an interruption. Use -j at most the number of physical cores.
Measured on an Apple M3 (one thread, 200 taxa, 1000 genes): wastral hybrid 851 s, tree-qmc hybrid 194 s.

Commands (GENES: the gene trees; OUT: the output tree; support options only with hybrid and support):
  wastral    BINARY -t 1 --mode 1|2|3|4 [-n MIN -x MAX -d MISSING] -o OUT GENES
             (--mode 1 hybrid, 2 support, 3 length, 4 none; ASTER's default search, -r 4 -s 4)
  tree-qmc   BINARY -i GENES -o OUT -w h|s|l|n [-n MIN -x MAX -d MISSING]
  wqfm-tree  BINARY -i GENES -o OUT -t 1 [--weight hybrid|support|length --support-scale MIN MAX]

Building the other methods (tested with ASTER ddc3dc6 and TREE-QMC v3.0.4):
  wastral    git clone https://github.com/chaoszhang/ASTER && cd ASTER && make wastral     (bin/wastral)
  tree-qmc   git clone --branch v3.0.4 https://github.com/molloy-lab/TREE-QMC && cd TREE-QMC/external/MQLib
             && make && cd ../.. && mkdir -p build && cd build && g++ -std=c++11 -O2
             -I ../external/MQLib/include -I ../external/toms743 -o tree-qmc ../src/*.cpp
             ../external/toms743/toms743.cpp ../external/MQLib/bin/MQLib.a -lm -DVERSION="\\"$(cat ../version.txt)\\""
             (TREE-QMC's main branch, 4.x, fails to build this way; its CMake build needs R)

Resuming (as in run_all.py): Without --overwrite, finished runs (species.tre and run.json present) are
skipped; interrupted runs start over. results.csv is repaired from finished runs.
With --overwrite, selected runs are executed again, including previously completed runs.
Previous selected results are removed before each rerun so failures are not mistaken for success.
A folder holds one method and setting: if settings.json differs from the options given, the script
stops (use another folder). Interrupting (Ctrl-C, or kill PID) stops the runs in progress. Two sweeps
never share a folder (a lock file guards it where the file system supports locks).

results.csv columns: method, weight, then those of run_all.py: condition, replicate, genes, taxa,
threads, seconds, max_rss_mb, fn_rate (true-tree branches missing from the estimate), fp_rate
(estimated branches not in the true tree), rf (normalised Robinson-Foulds distance), estimated_tree,
true_tree. Compare two folders with tools/weighted/compare_fn.py A/results.csv B/results.csv.
"""

import argparse
import concurrent.futures
import csv
import json
import os
import signal
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_all  # noqa: E402
from run_all import CONDITIONS, read_trees, write_atomic, compare, replicates, parse_list  # noqa: E402

METHODS = ('wastral', 'tree-qmc', 'wqfm-tree')
WEIGHTS = ('hybrid', 'support', 'length', 'none')
FIELDS = ['method', 'weight'] + run_all.FIELDS


def fmt(x):
    return f'{x:g}'


def command(method, weight, binary, genes, out, smin, smax, missing, extra):
    """The command line of one run."""
    supp = weight in ('hybrid', 'support')
    if method == 'wastral':
        mode = {'hybrid': '1', 'support': '2', 'length': '3', 'none': '4'}[weight]
        cmd = [binary, '-t', '1', '--mode', mode]
        if supp:
            cmd += ['-n', fmt(smin), '-x', fmt(smax), '-d', fmt(missing)]
        return cmd + extra + ['-o', out, genes]
    if method == 'tree-qmc':
        cmd = [binary, '-i', genes, '-o', out, '-w', {'hybrid': 'h', 'support': 's', 'length': 'l', 'none': 'n'}[weight]]
        if supp:
            cmd += ['-n', fmt(smin), '-x', fmt(smax), '-d', fmt(missing)]
        return cmd + extra
    cmd = [binary, '-i', genes, '-o', out, '-t', '1']
    if weight != 'none':
        cmd += ['--weight', weight]
    if supp:
        cmd += ['--support-scale', fmt(smin), fmt(smax)]
    return cmd + extra


class Runner:
    """Runs commands, each in its own process; stop() kills the ones still running."""

    def __init__(self, timeout):
        self.timeout = timeout
        self.procs = set()
        self.lock = threading.Lock()
        self.stopping = False
        self.started = {}                    # run folder -> start time, of the runs in progress

    def run(self, cmd, run_dir):
        """Run cmd in run_dir, logging to run.log; species.tre appears only if the run succeeds."""
        tree, tmp = os.path.join(run_dir, 'species.tre'), os.path.join(run_dir, 'species.tre.part')
        # a failed or interrupted rerun must never look completed because of old output
        for path in (tree, tmp, os.path.join(run_dir, 'run.json')):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
        env = dict(os.environ, OMP_NUM_THREADS='1')
        with open(os.path.join(run_dir, 'run.log'), 'w') as log:
            log.write('# ' + ' '.join(cmd) + '\n')
            log.flush()
            with self.lock:
                if self.stopping:
                    raise RuntimeError('stopped')
                t0 = time.perf_counter()
                p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                     cwd=run_dir, env=env)
                self.procs.add(p)
                self.started[run_dir] = time.time()
            try:
                timed_out = False
                while True:
                    pid, status, usage = os.wait4(p.pid, os.WNOHANG)
                    if pid:
                        break
                    if self.timeout and time.perf_counter() - t0 > self.timeout and not timed_out:
                        timed_out = True
                        p.kill()
                    time.sleep(0.25)
                secs = time.perf_counter() - t0
                p.returncode = os.waitstatus_to_exitcode(status)
            finally:
                with self.lock:
                    self.procs.discard(p)
                    self.started.pop(run_dir, None)
        if timed_out:
            raise RuntimeError(f'timed out after {self.timeout:g} s, see {run_dir}/run.log')
        if p.returncode != 0:
            raise RuntimeError(f'exited with {p.returncode}, see {run_dir}/run.log')
        trees = read_trees(tmp) if os.path.isfile(tmp) else []
        if not trees or not trees[0].endswith(';'):
            raise RuntimeError(f'no species tree written, see {run_dir}/run.log')
        rss = usage.ru_maxrss / (1 << 20) if sys.platform == 'darwin' else usage.ru_maxrss / 1024
        write_atomic(os.path.join(run_dir, 'run.json'),
                     json.dumps({'seconds': round(secs, 3), 'max_rss_mb': round(rss, 1), 'threads': 1}))
        os.replace(tmp, tree)
        return secs

    def stop(self):
        with self.lock:
            self.stopping = True
            for p in self.procs:
                if p.returncode is None:
                    try:
                        p.kill()
                    except OSError:
                        pass


def lock_folder(out):
    """Hold an exclusive lock on OUTPUT/.lock for the life of the process, where locks are supported."""
    try:
        import fcntl
    except ImportError:
        return None
    f = open(os.path.join(out, '.lock'), 'w')
    try:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit(f'another sweep is using {out}')
    except OSError:
        print(f'note: no file locking on this file system; do not start two sweeps on {out}', file=sys.stderr)
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('--method', required=True, choices=METHODS)
    ap.add_argument('--binary', required=True)
    ap.add_argument('--weight', default='hybrid', choices=WEIGHTS)
    ap.add_argument('-j', '--jobs', type=int, default=1)
    ap.add_argument('--genes')
    ap.add_argument('--conditions')
    ap.add_argument('--replicates')
    ap.add_argument('--gene-file', default='estimatedgenetre')
    ap.add_argument('--true-trees')
    ap.add_argument('--support-scale', nargs=2, type=float, default=[0.0, 1.0], metavar=('MIN', 'MAX'))
    ap.add_argument('--support-missing', type=float)
    ap.add_argument('--timeout', type=float, default=0)
    ap.add_argument('--extra', default='')
    ap.add_argument('--heartbeat', type=float, default=60)
    ap.add_argument('--overwrite', action='store_true')
    a = ap.parse_args()

    smin, smax = a.support_scale
    missing = smax if a.support_missing is None else a.support_missing
    if not smin < smax or not smin <= missing <= smax:
        sys.exit('need MIN < MAX and MIN <= --support-missing <= MAX')
    if a.method == 'wqfm-tree' and a.weight in ('hybrid', 'support') and missing != smax:
        sys.exit('wqfm-tree counts a branch without support as fully supported; --support-missing must be MAX')
    if a.jobs < 1:
        sys.exit('--jobs must be at least 1')
    binary = os.path.abspath(a.binary)
    if not (os.path.isfile(binary) and os.access(binary, os.X_OK)):
        sys.exit(f'no executable at {binary}')
    data, out = os.path.abspath(a.input), os.path.abspath(a.output)
    inp = os.path.join(data, run_all.GENE_TREES_DIR)
    true_root = os.path.abspath(a.true_trees) if a.true_trees else os.path.join(data, run_all.TRUE_TREES_DIR)
    if not os.path.isdir(inp):
        sys.exit(f'no {run_all.GENE_TREES_DIR} folder in {data}')
    if not os.path.isdir(true_root):
        print(f'warning: no true trees at {true_root}; accuracy columns stay empty', file=sys.stderr)
    os.makedirs(out, exist_ok=True)
    lock = lock_folder(out)  # noqa: F841 (held until exit)

    # one method and setting per folder
    supp = a.weight in ('hybrid', 'support')
    settings = {'method': a.method, 'weight': a.weight, 'extra': a.extra,
                'support_scale': [smin, smax] if supp else None, 'support_missing': missing if supp else None}
    settings_path = os.path.join(out, 'settings.json')
    if os.path.isfile(settings_path):
        old = json.load(open(settings_path))
        diff = [k for k in settings if old.get(k) != settings[k]]
        if diff:
            sys.exit(f'{out} holds runs with other settings ('
                     + ', '.join(f'{k}: {old.get(k)!r} there, {settings[k]!r} now' for k in diff)
                     + '); use another output folder')
        if old.get('binary') != binary:
            print(f'note: binary changed from {old.get("binary")} to {binary}', file=sys.stderr)
    write_atomic(settings_path, json.dumps(dict(settings, binary=binary), indent=1) + '\n')

    # every (condition, replicate, genes) to run, in order
    conds = a.conditions.split(',') if a.conditions else CONDITIONS
    for cond in conds:
        if cond not in CONDITIONS:
            sys.exit(f'unknown condition {cond}; expected one of:\n  ' + '\n  '.join(CONDITIONS))
    want_reps = {f'{r:02d}' for r in parse_list(a.replicates)} if a.replicates else None
    jobs = []
    for cond in conds:
        if not os.path.isdir(os.path.join(inp, cond)):
            print(f'skip {cond}: no folder {os.path.join(inp, cond)}', file=sys.stderr)
            continue
        for rep in replicates(os.path.join(inp, cond)):
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

    def complete(r):
        if None in r or any(r.get(k) is None for k in FIELDS):
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

    def row_for(cond, rep, genes, d):
        r = run_all.row_for(cond, rep, genes, d, true_tree(cond, rep))
        return dict(r, method=a.method, weight=a.weight)
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
        d = os.path.join(out, cond, rep, f'{genes or n}-genes')
        finished = os.path.isfile(os.path.join(d, 'species.tre')) and os.path.isfile(os.path.join(d, 'run.json'))
        if finished and not a.overwrite:
            if key not in rows:
                rows[key] = row_for(cond, rep, genes or n, d)
        else:
            # Remove old CSV data, so a forced rerun replaces rather than duplicates it.
            rows.pop(key, None)
            todo.append((cond, rep, genes, src, n, d))

    def save_csv():
        ordered = sorted(rows.values(), key=lambda r: (r['condition'], r['replicate'], int(r['genes'])))
        tmp = csv_path + '.part'
        with open(tmp, 'w', newline='') as f:
            w = csv.DictWriter(f, FIELDS, extrasaction='ignore')
            w.writeheader()
            w.writerows(ordered)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, csv_path)
    save_csv()
    done_before = len(jobs) - len(todo)
    print(f'{a.method} ({a.weight}): {len(jobs)} runs: {done_before} already finished, {len(todo)} to go; '
          f'{a.jobs} at a time, 1 thread each' + (' [overwrite enabled]' if a.overwrite else ''), flush=True)
    print('example command: ' + ' '.join(command(a.method, a.weight, binary, 'GENES', 'OUT', smin, smax, missing,
                                                 a.extra.split())), flush=True)

    runner = Runner(a.timeout)

    def stop_on_signal(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop_on_signal)

    started = [0]
    number = {}                             # run folder -> its number in the progress lines
    print_lock = threading.Lock()

    def job(cond, rep, genes, src, n, d):
        os.makedirs(d, exist_ok=True)
        gt = src
        if genes is not None and genes < n:               # the first `genes` trees of the file
            gt = os.path.join(d, 'gene-trees.tre')
            if a.overwrite or not os.path.isfile(gt):
                write_atomic(gt, '\n'.join(read_trees(src)[:genes]) + '\n')
        cmd = command(a.method, a.weight, binary, gt, os.path.join(d, 'species.tre.part'), smin, smax, missing,
                      a.extra.split())
        with print_lock:
            started[0] += 1
            number[d] = done_before + started[0]
            print(f'[{number[d]}/{len(jobs)}] starting {cond}/{rep}, {genes or n} genes '
                  f'({time.strftime("%H:%M:%S")})', flush=True)
        runner.run(cmd, d)
        return row_for(cond, rep, genes or n, d)

    failed = []
    t_start = time.time()
    finished = 0
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs)
    try:
        with open(csv_path, 'a', newline='') as f:
            w = csv.DictWriter(f, FIELDS, extrasaction='ignore')
            futures = {pool.submit(job, *t): t for t in todo}
            pending = set(futures)
            next_beat = time.time() + a.heartbeat
            while pending:
                wait = max(1.0, next_beat - time.time()) if a.heartbeat > 0 else None
                done, pending = concurrent.futures.wait(pending, timeout=wait,
                                                        return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    cond, rep, genes, src, n, d = futures[fut]
                    finished += 1
                    try:
                        row = fut.result()
                    except Exception as e:
                        failed.append(f'{cond}/{rep} ({genes or n} genes)')
                        with print_lock:
                            print(f'FAILED {cond}/{rep} ({genes or n} genes): {e}', file=sys.stderr, flush=True)
                        continue
                    rows[(cond, rep, str(genes or n))] = row
                    w.writerow(row)
                    f.flush()
                    os.fsync(f.fileno())
                    elapsed = time.time() - t_start
                    eta = elapsed / finished * (len(todo) - finished)
                    with print_lock:
                        print(f'[{number.get(d, "?")}/{len(jobs)}] done {cond}/{rep}, {genes or n} genes: '
                              f'{row["seconds"]} s, FN {row["fn_rate"] or "-"}, RF {row["rf"] or "-"} | '
                              f'elapsed {elapsed / 3600:.2f} h, estimated {eta / 3600:.2f} h left', flush=True)
                if a.heartbeat > 0 and time.time() >= next_beat:
                    next_beat = time.time() + a.heartbeat
                    with runner.lock:
                        running = sorted((t0, d) for d, t0 in runner.started.items())
                    with print_lock:
                        for t0, d in running:
                            print(f'    ... {time.time() - t0:6.0f} s  {os.path.relpath(d, out)}  '
                                  f'{run_all.last_line(os.path.join(d, "run.log"))}', flush=True)
    except KeyboardInterrupt:
        print('\ninterrupted: stopping the runs in progress (finished runs are kept; rerun to resume)',
              file=sys.stderr, flush=True)
        runner.stop()
        pool.shutdown(wait=True, cancel_futures=True)
        save_csv()
        sys.exit(130)
    pool.shutdown(wait=True)
    save_csv()                                               # sorted
    run_all.summary(rows.values(), os.path.join(out, 'summary.txt'))
    if failed:
        print(f'{len(failed)} run(s) failed (rerun to retry):\n  ' + '\n  '.join(failed), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
