import random, sys
sys.setrecursionlimit(100000)
import sim
n, k, H, seed, out = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
rnd = random.Random(seed)
nodes, root, h = sim.yule(n, rnd)
sp = sim.species_nested(nodes, root)
all_taxa = frozenset(str(i + 1) for i in range(n))
spb = sim.bipartitions(sp, all_taxa, '1')
rf = []
with open(out, 'w') as f:
    for g in range(k):
        t = sim.msc_gene_tree(nodes, root, h, H / h, rnd)
        f.write(sim.to_newick(t) + ';\n')
        if g < 20:
            gb = sim.bipartitions(t, all_taxa, '1'); rf.append(len(gb ^ spb) / (2 * (n - 3)))
print(f'{out}: n={n} k={k} H={H} mean RF(gene,species) over first 20 genes = {sum(rf)/len(rf):.3f}')
