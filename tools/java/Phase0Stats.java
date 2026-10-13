package src;

import java.util.ArrayList;
import java.util.Arrays;

import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.PreProcessing.GeneTrees;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * Phase-0 measurement (read-only; does not change the algorithm).
 * For a subproblem, measures how much the gene trees shrink if every maximal subtree whose
 * leaves all lie under a single dummy taxon is collapsed to one weighted leaf, and sibling
 * collapsed leaves of the same dummy at a node are merged. Purity is computed on the
 * unrooted tree (both below and above each node).
 */
public class Phase0Stats {

    static final int MIXED = -1, EMPTY = -2;

    static int combine(int a, int b) {
        if (a == EMPTY) return b;
        if (b == EMPTY) return a;
        return (a == b && a >= 0) ? a : MIXED;
    }

    public static String subproblem(GeneTrees gts, TaxaPerLevelWithPartition tx) {
        long fullNodes = 0, keptNodes = 0, fullDeg = 0, keptDeg = 0, fullCube = 0, keptCube = 0, fullSq = 0;
        long fullNodesRelevant = 0, polyFull = 0, polyKept = 0;
        int relevant = 0, unassigned = 0;
        int maxKeptDeg = 0, maxFullDeg = 0;

        for (Tree t : gts.geneTrees) {
            int N = t.nodes.size();
            int[] down = new int[N];
            int[] up = new int[N];
            boolean[] dummySeen = new boolean[tx.dummyTaxonCount];
            int units = 0;

            for (TreeNode v : t.topSortedNodes) {
                if (v.isLeaf()) {
                    int id = v.taxon.id;
                    if (tx.isInRealTaxa(id)) { down[v.index] = MIXED; units++; }
                    else if (tx.isInDummyTaxa(id)) {
                        int x = tx.inWhichDummyTaxa(id);
                        down[v.index] = x;
                        if (!dummySeen[x]) { dummySeen[x] = true; units++; }
                    } else { down[v.index] = EMPTY; unassigned++; }
                } else {
                    int acc = EMPTY;
                    for (TreeNode c : v.childs) acc = combine(acc, down[c.index]);
                    down[v.index] = acc;
                }
            }
            boolean isRelevant = units >= 4;
            if (isRelevant) relevant++;

            // top-down: up[child] = combine(up[parent], down[siblings])
            ArrayList<TreeNode> order = t.topSortedNodes;
            up[t.root.index] = EMPTY;
            for (int oi = order.size() - 1; oi >= 0; --oi) {
                TreeNode p = order.get(oi);
                if (p.isLeaf()) continue;
                int m = p.childs.size();
                int[] suffix = new int[m + 1];
                suffix[m] = EMPTY;
                for (int i = m - 1; i >= 0; --i) suffix[i] = combine(suffix[i + 1], down[p.childs.get(i).index]);
                int prefix = up[p.index];
                for (int i = 0; i < m; ++i) {
                    TreeNode c = p.childs.get(i);
                    up[c.index] = combine(prefix, suffix[i + 1]);
                    prefix = combine(prefix, down[c.index]);
                }
            }

            for (TreeNode u : order) {
                if (u.isLeaf() || u.isRoot()) continue;    // root is the artificial degree-2 node
                int m = u.childs.size();
                int D = m + 1;
                int[] labels = new int[D];
                for (int i = 0; i < m; ++i) labels[i] = down[u.childs.get(i).index];
                labels[m] = up[u.index];
                Arrays.sort(labels);
                int mixed = 0, distinctPure = 0, maxRun = 0, run = 0, prev = Integer.MIN_VALUE;
                for (int l : labels) {
                    if (l < 0) { mixed++; prev = Integer.MIN_VALUE; run = 0; continue; }
                    if (l == prev) run++; else { run = 1; distinctPure++; prev = l; }
                    maxRun = Math.max(maxRun, run);
                }
                fullNodes++;
                fullDeg += D;
                fullCube += (long) D * D * D;
                fullSq += (long) D * D;
                maxFullDeg = Math.max(maxFullDeg, D);
                if (D > 3) polyFull++;
                if (!isRelevant) continue;
                fullNodesRelevant++;
                boolean kept = maxRun < D - 1;
                if (!kept) continue;
                int Dp = mixed + distinctPure;
                keptNodes++;
                keptDeg += Dp;
                keptCube += (long) Dp * Dp * Dp;
                maxKeptDeg = Math.max(maxKeptDeg, Dp);
                if (Dp > 3) polyKept++;
            }
        }
        return String.format(
            "trees=%d relevant=%d fullNodes=%d keptNodes=%d fullDeg=%d keptDeg=%d fullCube=%d keptCube=%d fullSq=%d polyFull=%d polyKept=%d maxFullDeg=%d maxKeptDeg=%d unassigned=%d",
            gts.geneTrees.size(), relevant, fullNodes, keptNodes, fullDeg, keptDeg, fullCube, keptCube, fullSq,
            polyFull, polyKept, maxFullDeg, maxKeptDeg, unassigned);
    }
}
