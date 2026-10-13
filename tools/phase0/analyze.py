"""Aggregate Phase-0 instrumentation output (P0CAND / P0SUB lines) from an instrumented wQFM-TREE run.

Each subproblem's own P0CAND line is printed by its parent's divide (or by runWQFM for the root),
before the parent's P0SUB line; the recursion is depth-first, so a stack replays the pairing.

usage: analyze.py <run log> [n_taxa]   (n_taxa defaults to the "taxon count" line in the log)
"""
import math
import re
import sys
from collections import defaultdict

KV = re.compile(r'(\w+)=(-?\d+)')


def parse(path):
    subs, stack, pending, first, mismatches = [], [], [], True, 0
    for line in open(path):
        if line.startswith('P0CAND'):
            c = {k: int(v) for k, v in KV.findall(line)}
            if first:
                stack.append(c)
                first = False
            else:
                pending.append(c)
        elif line.startswith('P0SUB'):
            s = {k: int(v) for k, v in KV.findall(line)}
            own = stack.pop() if stack else None
            if own is None or (own['s'], own['r'], own['d']) != (s['s'], s['r'], s['d']):
                mismatches += 1
            s['cand'] = own
            subs.append(s)
            for c in reversed(pending):
                stack.append(c)
            pending = []
    return subs, mismatches


def main(path, n_taxa):
    subs, mism = parse(path)
    if not subs:
        print('no subproblems yet')
        return
    print(f'subproblems logged: {len(subs)}  (pairing mismatches: {mism})')

    fm_ms = sum(s['fmMs'] for s in subs)
    cand_ms = sum(s['cand']['ms'] for s in subs if s['cand'])
    print(f'measured time: FM passes {fm_ms/1000:.1f} s, consensus initial split {cand_ms/1000:.1f} s '
          f'({100*cand_ms/max(1, fm_ms+cand_ms):.0f}% of the two)')

    # FM work proxies: moves per pass ~ s; per move the Java touches every node of every gene tree
    jf = sum(s['passes'] * s['s'] * s['fullNodes'] for s in subs)
    cf = sum(s['passes'] * s['s'] * s['keptNodes'] for s in subs)
    # normalised so a binary node (D=3) costs 1 unit in both models: Java polytomy gains ~ (D/3)^3, closed form ~ D'/3
    jcube = sum(s['passes'] * s['s'] * s['fullCube'] for s in subs) / 27
    cdeg = sum(s['passes'] * s['s'] * s['keptDeg'] for s in subs) / 3
    # candidate scoring proxies: Java rescores with a full reset; proposal scores distinct splits on compressed trees
    jc = sum(s['cand']['rescored'] * s['fullNodes'] for s in subs if s['cand'])
    cc = sum(s['cand']['distinct'] * s['keptNodes'] for s in subs if s['cand'])
    cons = sum(s['cand']['considered'] for s in subs if s['cand'])
    resc = sum(s['cand']['rescored'] for s in subs if s['cand'])
    dist = sum(s['cand']['distinct'] for s in subs if s['cand'])

    print()
    print('FM phase (node visits; moves/pass taken as s):')
    print(f'  compression alone:                     {jf:.3e} -> {cf:.3e}  ({jf/max(1,cf):.1f}x fewer)')
    print(f'  + closed-form polytomies (binary=1):   {jcube:.3e} -> {cdeg:.3e}  ({jcube/max(1,cdeg):.1f}x fewer)')
    print('Initial-split phase:')
    print(f'  candidates considered {cons}, rescored by Java {resc}, distinct {dist}  '
          f'({resc/max(1,dist):.1f}x redundant rescoring)')
    print(f'  dedup + compression (node visits):    {jc:.3e} -> {cc:.3e}  ({jc/max(1,cc):.1f}x fewer)')

    print()
    print('By recursion level (weighted by s):')
    print(f'  {"level":>5} {"#sub":>5} {"mean s":>7} {"max d":>5} {"kept/full nodes":>16} {"relevant trees":>15} {"FM s":>8} {"cand s":>8}')
    lv = defaultdict(list)
    for s in subs:
        lv[s['level']].append(s)
    for L in sorted(lv):
        g = lv[L]
        ws = sum(x['s'] for x in g)
        kept = sum(x['s'] * x['keptNodes'] for x in g) / max(1, sum(x['s'] * x['fullNodes'] for x in g))
        rel = sum(x['relevant'] for x in g) / max(1, sum(x['trees'] for x in g))
        print(f'  {L:>5} {len(g):>5} {ws/len(g):>7.1f} {max(x["d"] for x in g):>5} {kept:>16.3f} {rel:>15.3f} '
              f'{sum(x["fmMs"] for x in g)/1000:>8.1f} {sum(x["cand"]["ms"] for x in g if x["cand"])/1000:>8.1f}')

    maxd = max(s['d'] for s in subs)
    bound = 2 * math.ceil(4 * math.log2(n_taxa))
    print()
    print(f'max dummies in a subproblem: {maxd}  (Java array bound before crash: {bound})')
    print(f'max polytomy degree: full {max(s["maxFullDeg"] for s in subs)}, compressed {max(s["maxKeptDeg"] for s in subs)}')
    print(f'FM passes per subproblem: mean {sum(s["passes"] for s in subs)/len(subs):.2f}, max {max(s["passes"] for s in subs)}')


def taxon_count(path):
    for line in open(path):
        if line.startswith('taxon count'):
            return int(line.split(':')[1])
    raise SystemExit('no "taxon count" line in log; pass n_taxa explicitly')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else taxon_count(sys.argv[1]))
