import java.io.*;
import java.util.*;

import src.Phase0Stats;
import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.PreProcessing.GeneTrees;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * Checks that collapsing maximal single-dummy subtrees (and merging same-dummy sibling clusters)
 * preserves Score(A,B,G) exactly, by brute-force quartet enumeration on both versions.
 * Also cross-checks Phase0Stats' kept-node count against an explicit construction.
 */
public class CompressCheck {

    static final int MIXED = -1, EMPTY = -2;

    static int combine(int a, int b) {
        if (a == EMPTY) return b;
        if (b == EMPTY) return a;
        return (a == b && a >= 0) ? a : MIXED;
    }

    // compressed leaf: unit = real taxon id (kind 0) or dummy index (kind 1)
    static class CLeaf { int node, kind, unit; double w; }

    static long keptTotal = 0;

    static double compressedScore(Tree t, TaxaPerLevelWithPartition tx) {
        int N = t.nodes.size();
        int[] down = new int[N], up = new int[N];
        double[] wdown = new double[N];
        Set<String> units = new HashSet<>();
        for (TreeNode v : t.topSortedNodes) {
            if (v.isLeaf()) {
                int id = v.taxon.id;
                wdown[v.index] = tx.getWeight(id);
                if (tx.isInRealTaxa(id)) { down[v.index] = MIXED; units.add("r" + id); }
                else { down[v.index] = tx.inWhichDummyTaxa(id); units.add("d" + tx.inWhichDummyTaxa(id)); }
            } else {
                int acc = EMPTY; double w = 0;
                for (TreeNode c : v.childs) { acc = combine(acc, down[c.index]); w += wdown[c.index]; }
                down[v.index] = acc; wdown[v.index] = w;
            }
        }
        double W = wdown[t.root.index];
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
        boolean relevant = units.size() >= 4;

        // kept nodes
        Map<TreeNode, Integer> id = new HashMap<>();
        for (TreeNode u : order) {
            if (u.isLeaf() || u.isRoot()) continue;
            int m = u.childs.size(), D = m + 1;
            Map<Integer, Integer> cnt = new HashMap<>();
            for (TreeNode c : u.childs) if (down[c.index] >= 0) cnt.merge(down[c.index], 1, Integer::sum);
            if (up[u.index] >= 0) cnt.merge(up[u.index], 1, Integer::sum);
            int maxRun = 0;
            for (int v : cnt.values()) maxRun = Math.max(maxRun, v);
            if (maxRun < D - 1) id.put(u, id.size());
        }
        if (relevant) keptTotal += id.size();
        if (id.isEmpty() || !relevant) return 0;

        List<List<Integer>> adj = new ArrayList<>();
        for (int i = 0; i < id.size(); i++) adj.add(new ArrayList<>());
        List<CLeaf> leaves = new ArrayList<>();
        Set<Long> edges = new HashSet<>();
        java.util.function.BiConsumer<Integer, Integer> addEdge = (a, b) -> {
            long key = ((long) Math.min(a, b) << 32) | Math.max(a, b);
            if (edges.add(key)) { adj.get(a).add(b); adj.get(b).add(a); }
        };
        java.util.function.BiFunction<Integer, CLeaf, Void> addLeaf = (u, lf) -> {
            int nid = adj.size(); adj.add(new ArrayList<>()); lf.node = nid; leaves.add(lf); addEdge.accept(u, nid); return null;
        };

        for (Map.Entry<TreeNode, Integer> e : id.entrySet()) {
            TreeNode u = e.getKey(); int uid = e.getValue();
            Map<Integer, Double> clusters = new TreeMap<>();
            for (TreeNode c : u.childs) {
                if (down[c.index] >= 0) clusters.merge(down[c.index], wdown[c.index], Double::sum);
                else if (c.isLeaf()) { CLeaf lf = new CLeaf(); lf.kind = 0; lf.unit = c.taxon.id; lf.w = 1; addLeaf.apply(uid, lf); }
                else { Integer cid = id.get(c); if (cid == null) throw new IllegalStateException("mixed child not kept"); addEdge.accept(uid, cid); }
            }
            if (up[u.index] >= 0) clusters.merge(up[u.index], W - wdown[u.index], Double::sum);
            else {
                TreeNode p = u.parent;
                TreeNode other = p;
                if (p.isRoot()) { other = p.childs.get(0) == u ? p.childs.get(1) : p.childs.get(0); }
                if (other.isLeaf()) {
                    if (!p.isRoot()) throw new IllegalStateException();
                    CLeaf lf = new CLeaf(); lf.kind = 0; lf.unit = other.taxon.id; lf.w = 1; addLeaf.apply(uid, lf);
                } else {
                    Integer oid = id.get(other); if (oid == null) throw new IllegalStateException("mixed parent not kept"); addEdge.accept(uid, oid);
                }
            }
            for (Map.Entry<Integer, Double> cl : clusters.entrySet()) {
                CLeaf lf = new CLeaf(); lf.kind = 1; lf.unit = cl.getKey(); lf.w = cl.getValue(); addLeaf.apply(uid, lf);
            }
        }

        // all-pairs leaf distances by BFS
        int L = leaves.size();
        int[][] dist = new int[L][];
        for (int i = 0; i < L; i++) {
            int[] d = new int[adj.size()]; Arrays.fill(d, -1);
            ArrayDeque<Integer> q = new ArrayDeque<>(); d[leaves.get(i).node] = 0; q.add(leaves.get(i).node);
            while (!q.isEmpty()) { int x = q.poll(); for (int y : adj.get(x)) if (d[y] < 0) { d[y] = d[x] + 1; q.add(y); } }
            dist[i] = new int[L];
            for (int j = 0; j < L; j++) dist[i][j] = d[leaves.get(j).node];
        }
        List<Integer> A = new ArrayList<>(), B = new ArrayList<>();
        for (int i = 0; i < L; i++) {
            CLeaf lf = leaves.get(i);
            int part = lf.kind == 0 ? tx.inWhichPartition(lf.unit) : tx.inWhichPartitionDummyTaxonByIndex(lf.unit);
            (part == 0 ? A : B).add(i);
        }
        double s = 0;
        for (int i = 0; i < A.size(); i++) for (int j = i + 1; j < A.size(); j++) {
            int a = A.get(i), b = A.get(j);
            if (same(leaves, a, b)) continue;
            for (int k = 0; k < B.size(); k++) for (int l = k + 1; l < B.size(); l++) {
                int c = B.get(k), d = B.get(l);
                if (same(leaves, c, d) || same(leaves, a, c) || same(leaves, a, d) || same(leaves, b, c) || same(leaves, b, d)) continue;
                double w = leaves.get(a).w * leaves.get(b).w * leaves.get(c).w * leaves.get(d).w;
                int s1 = dist[a][b] + dist[c][d], s2 = dist[a][c] + dist[b][d], s3 = dist[a][d] + dist[b][c];
                int tp = (s1 < s2 && s1 < s3) ? 1 : (s1 == s2 && s2 == s3) ? 0 : -1;
                s += w * tp;
            }
        }
        return s;
    }

    static boolean same(List<CLeaf> L, int a, int b) {
        return L.get(a).kind == 1 && L.get(b).kind == 1 && L.get(a).unit == L.get(b).unit;
    }

    public static void main(String[] args) throws Exception {
        long seed = Long.parseLong(args[0]);
        int trials = Integer.parseInt(args[1]);
        String tmpDir = args[2];
        BruteCheck.rnd = new Random(seed);
        Random rnd = BruteCheck.rnd;
        double maxRel = 0; int checks = 0, statMismatch = 0;
        long fullNodes = 0, kept = 0;
        for (int t = 0; t < trials; t++) {
            int n = 6 + rnd.nextInt(11);
            int k = 1 + rnd.nextInt(5);
            double collapse = rnd.nextDouble() < 0.5 ? 0 : 0.4;
            List<String> labels = new ArrayList<>();
            for (int i = 0; i < n; i++) labels.add("t" + i);
            File f = new File(tmpDir, "cc_" + seed + "_" + t + ".tre");
            try (PrintWriter pw = new PrintWriter(f)) {
                for (int g = 0; g < k; g++) {
                    List<String> present = new ArrayList<>(labels);
                    if (g > 0 && rnd.nextDouble() < 0.4) {
                        int drop = rnd.nextInt(Math.max(1, n - 5));
                        for (int i = 0; i < drop; i++) present.remove(rnd.nextInt(present.size()));
                    }
                    pw.println(BruteCheck.randomTree(present, collapse));
                }
            }
            RealTaxon.count = 0;
            GeneTrees gts = new GeneTrees(f.getPath());
            gts.readTaxaNames();
            gts.readGeneTrees(null);
            Object[] sp = BruteCheck.randomSubproblem(gts.taxa);
            RealTaxon[] rts = (RealTaxon[]) sp[0];
            DummyTaxon[] dts = (DummyTaxon[]) sp[1];
            if (rts.length + dts.length < 4) { f.delete(); continue; }
            TaxaPerLevelWithPartition tx = new TaxaPerLevelWithPartition(rts, dts, (int[]) sp[2], (int[]) sp[3], n);

            List<Map<Integer, int[]>> dists = new ArrayList<>();
            for (Tree tr : gts.geneTrees) dists.add(BruteCheck.leafDistances(tr, n));
            double full = BruteCheck.bruteScore(dists, tx, n);
            keptTotal = 0;
            double comp = 0;
            for (Tree tr : gts.geneTrees) comp += compressedScore(tr, tx);
            double rel = Math.abs(full - comp) / Math.max(1, Math.abs(full));
            maxRel = Math.max(maxRel, rel);
            if (rel > 1e-9) System.out.printf("MISMATCH trial=%d full=%.10f compressed=%.10f%n", t, full, comp);

            String st = Phase0Stats.subproblem(gts, tx);
            long statKept = Long.parseLong(st.replaceAll(".*keptNodes=(\\d+).*", "$1"));
            fullNodes += Long.parseLong(st.replaceAll(".*fullNodes=(\\d+).*", "$1"));
            kept += keptTotal;
            if (statKept != keptTotal) { statMismatch++; System.out.printf("STAT MISMATCH trial=%d stats=%d explicit=%d%n", t, statKept, keptTotal); }
            checks++;
            f.delete();
        }
        System.out.printf("seed=%d subproblems=%d  max rel err (full vs compressed)=%.2e  stat mismatches=%d  kept/full nodes=%d/%d%n",
            seed, checks, maxRel, statMismatch, kept, fullNodes);
        if (maxRel > 1e-9 || statMismatch > 0) System.exit(1);
    }
}
