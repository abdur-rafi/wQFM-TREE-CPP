"""Greedy (majority-rule extended) consensus of unrooted gene trees, as a PAUP-free stand-in.

usage: greedy_consensus.py <gene trees, one newick per line>  > consensus.tre
       greedy_consensus.py --clean <gene trees>               > topology-only trees (for the Java parser)

wQFM-TREE's run.sh uses PAUP 4.0b10 (a 32-bit Linux binary) for this step. Splits are added in
order of decreasing frequency (ties: first seen), keeping each one compatible with all accepted
splits. PAUP may break frequency ties differently, so initial bipartitions can differ in rare
cases; use PAUP's consensus when reproducing the paper exactly.
Only complete gene trees are supported (every tree must contain every taxon).
"""
import sys


def parse_newick(s):
    """Nested lists of leaf labels; ignores branch lengths, support values and [comments]."""
    pos = 0

    def skip_meta():
        nonlocal pos
        while pos < len(s) and s[pos] not in ',();':
            if s[pos] == '[':
                pos = s.index(']', pos) + 1
            else:
                pos += 1

    def node():
        nonlocal pos
        if s[pos] == '(':
            pos += 1
            kids = [node()]
            while s[pos] == ',':
                pos += 1
                kids.append(node())
            pos += 1                      # ')'
            skip_meta()                   # internal label / support / :length
            return kids
        start = pos
        while pos < len(s) and s[pos] not in ',():;[':
            pos += 1
        label = s[start:pos].strip().strip("'\"")
        skip_meta()
        return label

    s = s.strip()
    return node()


def leaves(t):
    return [t] if isinstance(t, str) else [x for c in t for x in leaves(c)]


def splits(t, all_taxa, ref):
    out = set()

    def rec(x):
        if isinstance(x, str):
            return frozenset([x])
        s = frozenset().union(*[rec(c) for c in x])
        side = s if ref not in s else all_taxa - s
        if 1 < len(side) < len(all_taxa) - 1:
            out.add(side)
        return s
    rec(t)
    return out


def greedy_consensus(trees):
    all_taxa = frozenset(leaves(trees[0]))
    for t in trees:
        if frozenset(leaves(t)) != all_taxa:
            raise SystemExit('incomplete gene tree found; this stand-in supports complete trees only')
    ref = min(all_taxa)
    counts, order = {}, []
    for t in trees:
        for s in splits(t, all_taxa, ref):
            if s not in counts:
                counts[s] = 0
                order.append(s)
            counts[s] += 1
    order.sort(key=lambda s: -counts[s])
    accepted = []
    for s in order:
        if all(s <= a or a <= s or not (s & a) for a in accepted):
            accepted.append(s)

    accepted.sort(key=len)
    children = {}
    for i, c in enumerate(accepted):
        parent = next((d for d in accepted[i + 1:] if c < d), None)
        children.setdefault(parent, []).append(c)

    def build(c):
        kids = children.get(c, [])
        covered = frozenset().union(*kids) if kids else frozenset()
        return [build(k) for k in kids] + sorted(c - covered)

    top = children.get(None, [])
    covered = frozenset().union(*top) if top else frozenset()
    return [build(k) for k in top] + sorted(all_taxa - covered - {ref}) + [ref]


def to_newick(t):
    return t if isinstance(t, str) else '(' + ','.join(to_newick(c) for c in t) + ')'


if __name__ == '__main__':
    clean = sys.argv[1] == '--clean'
    with open(sys.argv[2] if clean else sys.argv[1]) as f:
        trees = [parse_newick(line) for line in f if line.strip()]
    if clean:
        for t in trees:
            print(to_newick(t) + ';')
    else:
        print(to_newick(greedy_consensus(trees)) + ';')
