"""Simulate gene trees under the multispecies coalescent (see sim.py), optionally incomplete and with polytomies.

usage: gen_only.py <taxa> <genes> <height in coalescent units> <seed> <out gene trees>
                   [--missing P] [--collapse Q] [--species <out species tree>]

--missing P   drop each taxon from each gene tree with probability P (the first tree stays complete)
--collapse Q  contract each internal edge with probability Q (polytomies)
--species     also write the species tree (usable as the consensus tree for validation runs)
"""
import argparse
import random

import sim

ap = argparse.ArgumentParser()
for name, kind in [('n', int), ('k', int), ('H', float), ('seed', int), ('out', str)]:
    ap.add_argument(name, type=kind)
ap.add_argument('--missing', type=float, default=0.0)
ap.add_argument('--collapse', type=float, default=0.0)
ap.add_argument('--species')
a = ap.parse_args()

rnd = random.Random(a.seed)
nodes, root, h = sim.yule(a.n, rnd)
sp = sim.species_nested(nodes, root)
all_taxa = frozenset(str(i + 1) for i in range(a.n))
spb = sim.bipartitions(sp, all_taxa, '1')


def prune(t, keep):
    if isinstance(t, str):
        return t if t in keep else None
    kids = [c for c in (prune(c, keep) for c in t) if c is not None]
    if not kids:
        return None
    return kids[0] if len(kids) == 1 else kids


def collapse(t, q):
    if isinstance(t, str):
        return t
    out = []
    for c in t:
        c2 = collapse(c, q)
        if isinstance(c2, list) and rnd.random() < q:
            out.extend(c2)
        else:
            out.append(c2)
    return out


rf = []
with open(a.out, 'w') as f:
    for g in range(a.k):
        t = sim.msc_gene_tree(nodes, root, h, a.H / h, rnd)
        if g < 20:
            gb = sim.bipartitions(t, all_taxa, '1')
            rf.append(len(gb ^ spb) / (2 * (a.n - 3)))
        if a.missing and g > 0:
            t = prune(t, {x for x in all_taxa if rnd.random() >= a.missing})
        if a.collapse:
            t = collapse(t, a.collapse)
        f.write(sim.to_newick(t) + ';\n')
if a.species:
    with open(a.species, 'w') as f:
        f.write(sim.to_newick(sp) + ';\n')
print(f'{a.out}: n={a.n} k={a.k} H={a.H} missing={a.missing} collapse={a.collapse} '
      f'mean RF(gene,species) over first 20 genes = {sum(rf) / len(rf):.3f}')
