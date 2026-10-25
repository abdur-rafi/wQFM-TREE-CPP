"""Run a shared list of (method, weight, condition, replicate) jobs on several machines, resumably.

usage: python3 -I tools/astral2/run_queue.py --queue FILE --claims DIR -i DATA -o OUT [-j N] [--cores N]
                                              --bin wastral=PATH --bin astral=PATH --bin tree-qmc=PATH
                                              --bin wqfm-tree=PATH

  --queue FILE          One job per line, in order of priority (longest first packs best):
                          METHOD WEIGHT CONDITION REPLICATE [GENES] [genes=K] [threads=T]   ('#': comment)
                        METHOD and WEIGHT as in run_methods.py; GENES or genes=K: the first K gene trees
                        (default: all); threads=T: run with T threads (default 1; methods with -t only).
                        The file is read again whenever it changes, so jobs can be added or reordered
                        while the runners work.
  --claims DIR          Folder shared by all machines (e.g. on NFS). A machine claims a job by creating
                        DIR/<label>/<condition>/<replicate>-<genes> (mkdir is atomic; label: <method>-<weight>,
                        plus -t<T> for T > 1 threads), writes its host name there, and adds a file 'done' (or
                        'failed') when the run ends.
  -i, --input DIR       Data folder, as for run_all.py and run_methods.py
  -o, --output DIR      This machine's results: OUT/<label>/ is laid out as an output folder of run_methods.py
                        (results.csv, settings.json, <condition>/<replicate>/<genes>-genes/)
  -j, --jobs N          At most N runs of this runner at a time (default: --cores, or 1)
  --cores N             Cores of the machine, shared with other runners: a run of T threads takes min(T, N)
                        of them, and starts only when they are free. Busy cores are counted in /proc over all
                        processes of the queue's programs (their -t argument; 1 without one), so the runs of
                        other runners of this script count too. Whenever cores are free, the first unclaimed job
                        of the queue that fits starts; a job needing more cores waits until they are free.
  --min-threads T       Take only jobs with at least T threads (e.g. a machine kept for multi-threaded runs)
  --bin METHOD=PATH     Executable of a method (repeat for each method in the queue)
  --support-scale MIN MAX, --support-missing S, --timeout SECS, --heartbeat SECS: as in run_methods.py

Every machine runs this script with the same --queue and --claims and its own local --output; the folders of
all machines are merged afterwards (their runs are disjoint). Resuming: rerun the same command. Claims of this
machine without 'done' or 'failed' are released and run again; a finished run found in OUT is marked done
without rerunning. Remove a 'failed' claim folder to retry that job. Two runners on one machine must use
different --output folders.
"""

import argparse
import collections
import concurrent.futures
import csv
import json
import os
import signal
import socket
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_all  # noqa: E402
import run_methods  # noqa: E402
from run_all import read_trees, write_atomic  # noqa: E402

FIELDS = run_methods.FIELDS
Job = collections.namedtuple('Job', 'method weight cond rep genes src threads')


def label(job):
    return f'{job.method}-{job.weight}' + (f'-t{job.threads}' if job.threads > 1 else '')


def busy_cores(names, cores):
    """Cores taken by processes of the given programs on this machine that are not children of this
    process: each takes its -t argument (1 without one), at most `cores`."""
    me, n = os.getpid(), 0
    try:
        pids = [p for p in os.listdir('/proc') if p.isdigit()]
    except OSError:
        return 0
    for p in pids:
        try:
            stat = open(f'/proc/{p}/stat').read()
        except OSError:
            continue
        comm = stat[stat.index('(') + 1:stat.rindex(')')]
        ppid = int(stat[stat.rindex(')') + 2:].split()[1])
        if comm not in names or ppid == me:
            continue
        t = 1
        try:
            args = open(f'/proc/{p}/cmdline', 'rb').read().split(b'\0')
            for i, x in enumerate(args[:-1]):
                if x == b'-t':
                    t = int(args[i + 1])
        except (OSError, ValueError):
            pass
        n += min(max(t, 1), cores)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--queue', required=True)
    ap.add_argument('--claims', required=True)
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('-j', '--jobs', type=int)
    ap.add_argument('--cores', type=int)
    ap.add_argument('--min-threads', type=int, default=1)
    ap.add_argument('--bin', action='append', default=[], metavar='METHOD=PATH')
    ap.add_argument('--support-scale', nargs=2, type=float, default=[0.0, 1.0], metavar=('MIN', 'MAX'))
    ap.add_argument('--support-missing', type=float)
    ap.add_argument('--timeout', type=float, default=0)
    ap.add_argument('--heartbeat', type=float, default=600)
    a = ap.parse_args()
    max_jobs = a.jobs or a.cores or 1
    cores = a.cores or max_jobs

    smin, smax = a.support_scale
    missing = smax if a.support_missing is None else a.support_missing
    binaries = {}
    for b in a.bin:
        m, _, p = b.partition('=')
        if m not in run_methods.METHODS or not os.access(p, os.X_OK):
            sys.exit(f'bad --bin {b}: expected METHOD=EXECUTABLE with METHOD in {run_methods.METHODS}')
        binaries[m] = os.path.abspath(p)
    names = {os.path.basename(p)[:15] for p in binaries.values()}
    data, out, claims = os.path.abspath(a.input), os.path.abspath(a.output), os.path.abspath(a.claims)
    gene_root = os.path.join(data, run_all.GENE_TREES_DIR)
    true_root = os.path.join(data, run_all.TRUE_TREES_DIR)
    host = socket.gethostname()

    gene_counts, true_cache = {}, {}

    def gene_count(src):
        if src not in gene_counts:
            gene_counts[src] = len(read_trees(src))
        return gene_counts[src]

    def true_tree(cond, rep):
        if (cond, rep) not in true_cache:
            p = os.path.join(true_root, cond, rep, 's_tree.trees')
            true_cache[(cond, rep)] = read_trees(p)[0] if os.path.isfile(p) else None
        return true_cache[(cond, rep)]

    def paths(job):
        g = job.genes or gene_count(job.src)
        return (os.path.join(claims, label(job), job.cond, f'{job.rep}-{g}'),
                os.path.join(out, label(job), job.cond, job.rep, f'{g}-genes'), g)

    def finished_run(d):
        return os.path.isfile(os.path.join(d, 'species.tre')) and os.path.isfile(os.path.join(d, 'run.json'))

    def read_queue():
        jobs = []
        for ln, line in enumerate(open(a.queue), 1):
            f = line.split('#')[0].split()
            if not f:
                continue
            where = f'{a.queue}:{ln}'
            if len(f) < 4 or f[0] not in run_methods.METHODS or f[1] not in run_methods.WEIGHTS:
                raise ValueError(f'{where}: expected METHOD WEIGHT CONDITION REPLICATE [GENES] [genes=K] [threads=T]')
            genes, threads = None, 1
            for x in f[4:]:
                if x.isdigit():
                    genes = int(x)
                elif x.startswith('genes='):
                    genes = int(x[6:])
                elif x.startswith('threads='):
                    threads = int(x[8:])
                else:
                    raise ValueError(f'{where}: unknown field {x}')
            if f[0] not in binaries:
                raise ValueError(f'{where}: no --bin for {f[0]}')
            if threads < 1 or (threads > 1 and f[0] not in run_methods.THREADED):
                raise ValueError(f'{where}: {f[0]} cannot run with {threads} threads')
            if f[0] == 'astral' and f[1] != 'none':
                raise ValueError(f'{where}: astral is unweighted (weight none)')
            src = os.path.join(gene_root, f[2], f[3], 'estimatedgenetre')
            if not os.path.isfile(src):
                raise ValueError(f'{where}: no gene trees at {src}')
            jobs.append(Job(f[0], f[1], f[2], f[3], genes, src, threads))
        return jobs

    # per label: settings.json, and results.csv rebuilt from the finished runs of this machine
    csv_rows, files, writers = {}, {}, {}
    supp = lambda w: w in ('hybrid', 'support')

    def ensure_label(lab, job0, jobs):
        if lab in csv_rows:
            return
        mdir = os.path.join(out, lab)
        os.makedirs(mdir, exist_ok=True)
        w = job0.weight
        write_atomic(os.path.join(mdir, 'settings.json'), json.dumps(
            {'method': job0.method, 'weight': w, 'threads': job0.threads, 'extra': '',
             'support_scale': [smin, smax] if supp(w) else None, 'support_missing': missing if supp(w) else None,
             'binary': binaries[job0.method]}, indent=1) + '\n')
        rows = {}
        for job in jobs:
            if label(job) == lab:
                claim, d, g = paths(job)
                if finished_run(d):
                    r = run_all.row_for(job.cond, job.rep, g, d, true_tree(job.cond, job.rep))
                    rows[(job.cond, job.rep, str(g))] = dict(r, method=job.method, weight=w)
        path = os.path.join(mdir, 'results.csv')
        tmp = path + '.part'
        with open(tmp, 'w', newline='') as f:
            wr = csv.DictWriter(f, FIELDS, extrasaction='ignore')
            wr.writeheader()
            wr.writerows(sorted(rows.values(), key=lambda r: (r['condition'], r['replicate'], int(r['genes']))))
        os.replace(tmp, path)
        csv_rows[lab] = rows
        files[lab] = open(path, 'a', newline='')
        writers[lab] = csv.DictWriter(files[lab], FIELDS, extrasaction='ignore')

    def ensure_labels(jobs):
        first = {}
        for job in jobs:
            first.setdefault(label(job), job)
        for lab in sorted(first):
            ensure_label(lab, first[lab], jobs)

    try:
        jobs = read_queue()
    except ValueError as e:
        sys.exit(str(e))
    queue_mtime = os.path.getmtime(a.queue)
    ensure_labels(jobs)

    # release this machine's unfinished claims; mark runs finished here but not yet marked
    for job in jobs:
        claim, d, g = paths(job)
        if not os.path.isdir(claim):
            continue
        try:
            owner = open(os.path.join(claim, 'host')).read().strip()
        except OSError:
            owner = ''
        ended = os.path.exists(os.path.join(claim, 'done')) or os.path.exists(os.path.join(claim, 'failed'))
        if owner == host and not ended:
            if finished_run(d):
                open(os.path.join(claim, 'done'), 'w').close()
            else:
                for name in os.listdir(claim):
                    os.remove(os.path.join(claim, name))
                os.rmdir(claim)

    known = set()                             # claims (by any machine) seen: never tried again
    slots = lambda job: min(job.threads, cores)

    def next_job(free):
        """Claim the first unclaimed job of the queue that needs at most `free` cores. Returns (job, position),
        or None, with `pending` set to whether unclaimed jobs remain (that need more cores)."""
        nonlocal jobs, queue_mtime, pending
        try:
            mt = os.path.getmtime(a.queue)
            if mt != queue_mtime:
                new = read_queue()
                ensure_labels(new)
                jobs, queue_mtime = new, mt
                print(f'{time.strftime("%H:%M:%S")} queue reloaded: {len(jobs)} jobs', flush=True)
        except (OSError, ValueError) as e:
            print(f'{time.strftime("%H:%M:%S")} queue not reloaded: {e}', file=sys.stderr, flush=True)
        pending = False
        for pos, job in enumerate(jobs, 1):
            claim, d, g = paths(job)
            if claim in known or job.threads < a.min_threads:
                continue
            if os.path.isdir(claim):          # claimed elsewhere
                known.add(claim)
                continue
            if slots(job) > free:
                pending = True
                continue
            os.makedirs(os.path.dirname(claim), exist_ok=True)
            try:
                os.mkdir(claim)
            except FileExistsError:
                known.add(claim)
                continue
            known.add(claim)
            with open(os.path.join(claim, 'host'), 'w') as f:
                f.write(host + '\n')
            return job, pos
        return None

    runner = run_methods.Runner(a.timeout)
    lock = threading.Lock()

    def run(job):
        claim, d, g = paths(job)
        os.makedirs(d, exist_ok=True)
        gt = job.src
        if job.genes is not None and job.genes < gene_count(job.src):
            gt = os.path.join(d, 'gene-trees.tre')
            write_atomic(gt, '\n'.join(read_trees(job.src)[:job.genes]) + '\n')
        if not finished_run(d):
            cmd = run_methods.command(job.method, job.weight, binaries[job.method], gt,
                                      os.path.join(d, 'species.tre.part'), smin, smax, missing, [], job.threads)
            runner.run(cmd, d, job.threads)
        r = run_all.row_for(job.cond, job.rep, g, d, true_tree(job.cond, job.rep))
        return dict(r, method=job.method, weight=job.weight)

    def free_cores():
        mine = sum(slots(job) for job, t0 in running.values())
        n = cores - mine
        if a.cores:
            n -= busy_cores(names, cores)
        return n if len(running) < max_jobs else 0

    def stop_on_signal(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop_on_signal)

    print(f'{host}: {len(jobs)} jobs in the queue; up to {max_jobs} at a time on {cores} cores'
          + (' shared with other runners' if a.cores else ''), flush=True)
    t_start = time.time()
    done_here = failed_here = 0
    running = {}
    pending = False
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=max_jobs)
    next_beat = time.time() + a.heartbeat
    try:
        while True:
            exhausted = False
            while True:
                free = free_cores()
                if free <= 0:
                    break
                got = next_job(free)
                if got is None:
                    exhausted = not pending
                    break
                job, pos = got
                running[pool.submit(run, job)] = (job, time.time())
                print(f'{time.strftime("%H:%M:%S")} start {label(job)} {job.cond}/{job.rep} '
                      f'(queue position {pos}/{len(jobs)})', flush=True)
            if not running and exhausted:
                break
            if not running:                   # waiting for cores taken by other runners
                time.sleep(30)
                continue
            done, _ = concurrent.futures.wait(list(running), timeout=30,
                                              return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in done:
                job, t0 = running.pop(fut)
                claim, d, g = paths(job)
                lab = label(job)
                try:
                    row = fut.result()
                except Exception as e:
                    failed_here += 1
                    with open(os.path.join(claim, 'failed'), 'w') as f:
                        f.write(f'{e}\n')
                    print(f'{time.strftime("%H:%M:%S")} FAILED {lab} {job.cond}/{job.rep}: {e}', file=sys.stderr, flush=True)
                    continue
                with lock:
                    csv_rows[lab][(job.cond, job.rep, str(g))] = row
                    writers[lab].writerow(row)
                    files[lab].flush()
                    os.fsync(files[lab].fileno())
                open(os.path.join(claim, 'done'), 'w').close()
                done_here += 1
                print(f'{time.strftime("%H:%M:%S")} done {lab} {job.cond}/{job.rep}: {row["seconds"]} s, '
                      f'FN {row["fn_rate"] or "-"} | {done_here} done here, {len(running)} running, '
                      f'{(time.time() - t_start) / 3600:.2f} h elapsed', flush=True)
            if a.heartbeat > 0 and time.time() >= next_beat:
                next_beat = time.time() + a.heartbeat
                for job, t0 in sorted(running.values(), key=lambda x: x[1]):
                    claim, d, g = paths(job)
                    print(f'    ... {time.time() - t0:7.0f} s  {label(job)} {job.cond}/{job.rep}  '
                          f'{run_all.last_line(os.path.join(d, "run.log"))[:100]}', flush=True)
    except KeyboardInterrupt:
        print('\ninterrupted: stopping the runs in progress; rerun the same command to resume',
              file=sys.stderr, flush=True)
        runner.stop()
        pool.shutdown(wait=True, cancel_futures=True)
        sys.exit(130)
    pool.shutdown(wait=True)
    for lab, rows in csv_rows.items():
        if rows:
            print(f'== {lab} (this machine)')
            run_all.summary(rows.values(), os.path.join(out, lab, 'summary.txt'))
    print(f'{host}: queue exhausted; {done_here} runs done here, {failed_here} failed', flush=True)


if __name__ == '__main__':
    main()
