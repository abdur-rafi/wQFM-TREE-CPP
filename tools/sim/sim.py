"""Synthetic stand-in for ASTRAL-II-style data (no external dependencies).

Yule species tree scaled to height H (coalescent units) -> k gene trees under the multispecies
coalescent -> optional estimation-error proxy (random NNIs) -> greedy (majority-rule extended)
consensus. Writes gene trees, species tree, consensus, and prints mean normalised RF
(gene tree vs species tree) so the discordance level can be compared with ASTRAL-II conditions.
"""
import random
import sys


def yule(n, rnd):
    # nodes: dict id -> (children list, time); leaves time 0; build backwards from present
    lineages = [('L', i + 1) for i in range(n)]
    t = 0.0
    nodes = {}
    nid = 0
    # coalescent-like backwards Yule: with j lineages, waiting time ~ Exp(j) (rate-1 Yule)
    active = [(('leaf', i + 1), 0.0) for i in range(n)]
    while len(active) > 1:
        j = len(active)
        t += rnd.expovariate(j)
        a, b = rnd.sample(range(j), 2)
        x, y = active[a], active[b]
        nid += 1
        node = ('int', nid, [x[0], y[0]], t)
        nodes[nid] = node
        active = [active[i] for i in range(j) if i not in (a, b)] + [(('int', nid), t)]
    return nodes, active[0][0], t


def msc_gene_tree(nodes, root, height, scale, rnd):
    """Return nested-list gene tree. Species node times are multiplied by scale (CU)."""
    def node_time(ref):
        return 0.0 if ref[0] == 'leaf' else nodes[ref[1]][3] * scale

    def lineages_at_top(ref, parent_time):
        # lineages entering this species branch from below, then coalesce up to parent_time
        if ref[0] == 'leaf':
            lins = [str(ref[1])]
        else:
            _, nid, kids, _ = nodes[ref[1]]
            lins = []
            for k in kids:
                lins.extend(lineages_at_top(k, node_time(ref)))
        t = node_time(ref)
        while len(lins) > 1:
            j = len(lins)
            w = rnd.expovariate(j * (j - 1) / 2)
            if parent_time is not None and t + w > parent_time:
                break
            t += w
            a, b = rnd.sample(range(j), 2)
            merged = [lins[a], lins[b]]
            lins = [lins[i] for i in range(j) if i not in (a, b)] + [merged]
        return lins

    top = lineages_at_top(root, None)
    assert len(top) == 1
    return top[0]


def to_newick(t):
    return t if isinstance(t, str) else '(' + ','.join(to_newick(c) for c in t) + ')'


def species_nested(nodes, ref):
    if ref[0] == 'leaf':
        return str(ref[1])
    return [species_nested(nodes, k) for k in nodes[ref[1]][2]]


def leaves(t):
    return [t] if isinstance(t, str) else [x for c in t for x in leaves(c)]


def bipartitions(t, all_taxa, ref_taxon):
    """Set of nontrivial splits, each as frozenset side NOT containing ref_taxon."""
    out = set()

    def rec(x):
        if isinstance(x, str):
            return frozenset([x])
        s = frozenset().union(*[rec(c) for c in x])
        side = s if ref_taxon not in s else all_taxa - s
        if 1 < len(side) < len(all_taxa) - 1:
            out.add(side)
        return s
    rec(t)
    return out


def nni(t, rnd):
    """One random NNI on a rooted nested-list binary tree (treating it as unrooted)."""
    internal = []

    def collect(x, parent, idx):
        if isinstance(x, list):
            if parent is not None:
                internal.append((parent, idx))
            for i, c in enumerate(x):
                collect(c, x, i)
    collect(t, None, None)
    if not internal:
        return
    parent, idx = rnd.choice(internal)
    u = parent[idx]
    sib_idx = 1 - idx if len(parent) == 2 else rnd.choice([i for i in range(len(parent)) if i != idx])
    ci = rnd.randrange(len(u))
    u[ci], parent[sib_idx] = parent[sib_idx], u[ci]


def greedy_consensus(trees, all_taxa, ref):
    counts = {}
    order = []
    for t in trees:
        for s in bipartitions(t, all_taxa, ref):
            if s not in counts:
                counts[s] = 0
                order.append(s)
            counts[s] += 1
    order.sort(key=lambda s: -counts[s])          # stable: ties keep first-seen order
    accepted = []
    for s in order:
        if all(s <= a or a <= s or not (s & a) for a in accepted):
            accepted.append(s)
    # build tree rooted at ref: clades = accepted sides (none contain ref)
    accepted.sort(key=len)
    children_of = {}
    clade_nodes = {}
    for c in accepted:
        clade_nodes[c] = []
    taxa_parent = {}
    for i, c in enumerate(accepted):
        parent = None
        for d in accepted[i + 1:]:
            if c < d:
                parent = d
                break
        children_of.setdefault(parent, []).append(c)

    def build(c):
        kids = children_of.get(c, [])
        covered = frozenset().union(*kids) if kids else frozenset()
        parts = [build(k) for k in kids] + sorted(c - covered)
        return parts

    top_kids = children_of.get(None, [])
    covered = frozenset().union(*top_kids) if top_kids else frozenset()
    root = [build(k) for k in top_kids] + sorted((all_taxa - covered) - {ref}) + [ref]
    return root


def main():
    n, k, H, nni_frac, seed, out_prefix = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    rnd = random.Random(seed)
    nodes, root, h = yule(n, rnd)
    scale = H / h
    sp = species_nested(nodes, root)
    all_taxa = frozenset(str(i + 1) for i in range(n))
    sp_bip = bipartitions(sp, all_taxa, '1')
    gts = []
    rf = 0.0
    for _ in range(k):
        g = msc_gene_tree(nodes, root, h, scale, rnd)
        for _ in range(int(nni_frac * n)):
            nni(g, rnd)
        gts.append(g)
        gb = bipartitions(g, all_taxa, '1')
        rf += len(gb ^ sp_bip) / (2 * (n - 3))
    with open(out_prefix + '.gt.tre', 'w') as f:
        for g in gts:
            f.write(to_newick(g) + ';\n')
    with open(out_prefix + '.sp.tre', 'w') as f:
        f.write(to_newick(sp) + ';\n')
    cons = greedy_consensus(gts, all_taxa, '1')
    with open(out_prefix + '.cons.tre', 'w') as f:
        f.write(to_newick(cons) + ';\n')
    cb = bipartitions(cons, all_taxa, '1')
    print(f'{out_prefix}: n={n} k={k} H={H} nni/n={nni_frac} mean RF(gene,species)={rf/k:.3f} '
          f'consensus resolved={len(cb)}/{n-3} RF(cons,species)={len(cb ^ sp_bip)/(2*(n-3)):.3f}')


if __name__ == '__main__':
    main()
