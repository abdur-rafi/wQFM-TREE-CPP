"""Compare a Java trace (build/java-trace, -Dwqfm.trace=...) with a C++ trace (wqfm-tree --trace).

usage: trace_compare.py <java trace> <c++ trace> [tolerance, default 1e-9]

C++ replay traces (wqfm-tree --replay <java trace>) follow Java's decisions, so every candidate score
and every gain can be compared; where the C++ code would have chosen a different move, the gap between
the two gains tells a near-tie (rounding) from a real disagreement.
Normal traces are compared subproblem by subproblem until they diverge.
Differences are relative to the largest magnitude in the candidate set or FM pass (at least 1): rounding
error scales with it, and a value near zero carries the error of the large terms that cancel in it.
Exits non-zero on a disagreement beyond the tolerance or a structural mismatch.
"""
import sys


def parse(path):
    parts, subs, order, taxa = {}, {}, [], []
    cur_p = cur_s = None
    open_pass = False
    for line in open(path):
        line = line.rstrip('\n')
        if len(line) < 2 or line[1] != ' ':
            continue
        t, rest = line[0], line[2:]
        if t == 'T':
            taxa.append(rest)
        elif t == 'P':
            cur_p = parts[rest] = {'cand': [], 'R': None}
            order.append(('P', rest))
        elif t == 'C':
            cur_p['cand'].append(float(rest))
        elif t == 'R':
            cur_p['R'] = rest
        elif t == 'S':
            cur_s = subs[rest] = {'I': None, 'passes': [], 'F': None}
            open_pass = False
            order.append(('S', rest))
        elif t == 'I':
            cur_s['I'] = rest
        elif t == 'M':
            f = rest.split()
            if not open_pass:
                cur_s['passes'].append({'moves': [], 'E': None})
                open_pass = True
            mv = {'d': int(f[0]), 'i': int(f[1]), 'gain': float(f[2]), 'cg': float(f[3])}
            if len(f) >= 7:
                mv['own'] = (int(f[4]), int(f[5]), float(f[6]))
            cur_s['passes'][-1]['moves'].append(mv)
        elif t == 'E':
            if not open_pass:
                cur_s['passes'].append({'moves': [], 'E': None})
            cur_s['passes'][-1]['E'] = int(rest)
            open_pass = False
        elif t == 'F':
            cur_s['F'] = rest
    return parts, subs, order, taxa


def same_split(x, y):
    """bits equal up to swapping the two sides"""
    if x is None or y is None:
        return x == y
    return x == y or ''.join({'0': '1', '1': '0'}.get(c, c) for c in x) == y


def rel(a, b, scale):
    return abs(a - b) / max(1.0, abs(b), scale)


def main(jpath, cpath, tol):
    jparts, jsubs, jorder, jtaxa = parse(jpath)
    cparts, csubs, corder, ctaxa = parse(cpath)
    replay = any('own' in m for s in csubs.values() for p in s['passes'] for m in p['moves'])
    problems = []
    if jtaxa != ctaxa:
        problems.append('taxon order differs (pass --taxa-order <java trace> to the C++ run)')

    n_cand = n_moves = matched_s = matched_p = 0
    worst_cand = worst_gain = 0.0
    ties = disagreements = 0
    worst_disagree = 0.0
    first_divergence = None

    for kind, key in jorder:
        if kind == 'P':
            if key not in cparts:
                first_divergence = first_divergence or ('initial split of', key)
                continue
            matched_p += 1
            jc, cc = jparts[key]['cand'], cparts[key]['cand']
            if len(jc) != len(cc):
                problems.append(f'candidate count {len(jc)} vs {len(cc)} for P {key[:80]}')
            scale = max((abs(x) for x in jc), default=0.0)
            for a, b in zip(jc, cc):
                n_cand += 1
                worst_cand = max(worst_cand, rel(b, a, scale))
            if not same_split(jparts[key]['R'], cparts[key]['R']):
                first_divergence = first_divergence or ('initial split chosen for', key)
            continue

        if key not in csubs:
            first_divergence = first_divergence or ('subproblem', key)
            continue
        matched_s += 1
        js, cs = jsubs[key], csubs[key]
        if not same_split(js['I'], cs['I']):
            first_divergence = first_divergence or ('initial split of', key)
            continue
        for pi, (jp, cp) in enumerate(zip(js['passes'], cs['passes'])):
            scale = max((abs(m['gain']) for m in jp['moves']), default=0.0)
            for mi, (jm, cm) in enumerate(zip(jp['moves'], cp['moves'])):
                if (jm['d'], jm['i']) != (cm['d'], cm['i']):
                    first_divergence = first_divergence or (f'FM pass {pi} move {mi} (Java {jm["d"]}:{jm["i"]} gain {jm["gain"]:.6g}, '
                                                            f'C++ {cm["d"]}:{cm["i"]} gain {cm["gain"]:.6g}) in', key)
                    break
                n_moves += 1
                worst_gain = max(worst_gain, rel(cm['gain'], jm['gain'], scale))
                if 'own' in cm and cm['own'][:2] != (cm['d'], cm['i']):
                    gap = rel(cm['own'][2], cm['gain'], scale)
                    if gap <= tol:
                        ties += 1
                    else:
                        disagreements += 1
                        worst_disagree = max(worst_disagree, gap)
            if len(jp['moves']) != len(cp['moves']) or jp['E'] != cp['E']:
                first_divergence = first_divergence or (f'end of FM pass {pi} in', key)
        if len(js['passes']) != len(cs['passes']):
            first_divergence = first_divergence or ('number of FM passes in', key)
        if not same_split(js['F'], cs['F']):
            first_divergence = first_divergence or ('final split of', key)

    print(f'mode: {"replay (C++ follows Java decisions)" if replay else "independent runs"}')
    print(f'subproblems matched: {matched_s}/{sum(1 for k, _ in jorder if k == "S")}, '
          f'initial-split computations matched: {matched_p}/{sum(1 for k, _ in jorder if k == "P")}')
    print(f'candidate scores compared: {n_cand}, max rel diff {worst_cand:.2e}')
    print(f'FM gains compared: {n_moves}, max rel diff {worst_gain:.2e}')
    if replay:
        print(f'C++ would choose another move: {ties} near-ties (rel gap <= {tol:g}), '
              f'{disagreements} real disagreements (max rel gap {worst_disagree:.2e})')
    if first_divergence:
        what, key = first_divergence
        print(f'first divergence: {what} {key[:120]}')
    for p in problems:
        print('PROBLEM:', p)

    bad = problems or worst_cand > tol or worst_gain > tol or disagreements or (replay and first_divergence)
    print('RESULT:', 'MISMATCH' if bad else ('agreement' if not first_divergence else 'agreement up to the first divergence'))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 1e-9)
