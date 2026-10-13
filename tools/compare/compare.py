"""Exact comparison of Java (TopLevelBench) and C++ (gainsweep) outputs written under one prefix.

usage: compare.py <prefix>   (reads <prefix>.{java,cpp}.{gains,fmlog,fmpart})
Exits non-zero unless gains, FM moves and final partitions are identical.
"""
import sys


def table(path, conv):
    with open(path) as f:
        return [tuple(conv(x) for x in line.split()) for line in f if line.strip()]


def main(p):
    ok = True
    jg = dict(table(p + '.java.gains', float))
    cg = dict(table(p + '.cpp.gains', float))
    same = sum(1 for t in jg if cg.get(t) == jg[t])
    ok &= same == len(jg) == len(cg)
    print(f'gains identical:            {same}/{len(jg)}')

    jl, cl = table(p + '.java.fmlog', float), table(p + '.cpp.fmlog', float)
    same = sum(1 for a, b in zip(jl, cl) if a == b)
    ok &= same == len(jl) == len(cl)
    print(f'FM moves identical:         {same}/{max(len(jl), len(cl))}  (taxon, gain, cumulative gain)')

    jp, cp = table(p + '.java.fmpart', int), table(p + '.cpp.fmpart', int)
    ok &= jp == cp
    print(f'final partition identical:  {jp == cp}')
    print('RESULT:', 'exact agreement' if ok else 'MISMATCH')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main(sys.argv[1])
