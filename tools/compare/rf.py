"""Robinson-Foulds distance between two unrooted trees (first tree of each file).

usage: rf.py <tree file 1> <tree file 2>
Prints the number of splits only in each tree and the normalised RF distance; exits 0 iff identical.
"""
import sys


def parse(s):
    s = s.strip()
    pos = 0

    def node():
        nonlocal pos
        if s[pos] == '(':
            pos += 1
            kids = [node()]
            while s[pos] == ',':
                pos += 1
                kids.append(node())
            pos += 1
            while pos < len(s) and s[pos] not in ',();':
                pos += 1
            return kids
        start = pos
        while pos < len(s) and s[pos] not in ',();':      # label and optional :length
            pos += 1
        label = s[start:pos].split(':')[0].strip()
        return label

    return node()


def splits(t):
    taxa = set()

    def leaves(x):
        if isinstance(x, str):
            taxa.add(x)
            return frozenset([x])
        return frozenset().union(*map(leaves, x))
    leaves(t)
    ref = min(taxa)
    out = set()

    def rec(x):
        if isinstance(x, str):
            return frozenset([x])
        s = frozenset().union(*map(rec, x))
        side = s if ref not in s else frozenset(taxa) - s
        if 1 < len(side) < len(taxa) - 1:
            out.add(side)
        return s
    rec(t)
    return out, taxa


def first_tree(path):
    with open(path) as f:
        return next(line for line in f if '(' in line)


if __name__ == '__main__':
    (s1, t1), (s2, t2) = splits(parse(first_tree(sys.argv[1]))), splits(parse(first_tree(sys.argv[2])))
    if t1 != t2:
        sys.exit('trees have different taxon sets')
    n = len(t1)
    print(f'splits only in first: {len(s1 - s2)}, only in second: {len(s2 - s1)}, '
          f'normalised RF: {len(s1 ^ s2) / (2 * (n - 3)):.4f}')
    sys.exit(0 if s1 == s2 else 1)
