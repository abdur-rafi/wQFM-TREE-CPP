import java.io.*;
import java.util.*;

import src.DSPerLevel.BookKeepingPerLevelv2;
import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.PreProcessing.GeneTrees;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * Brute-force verification of wQFM-TREE's score and gain computation.
 * Score(A,B,G) = sum_g sum_{q: 2 taxa in FA, 2 in FB, no two under same dummy} w(q) * (+1 sat, -1 viol, 0 unresolved)
 */
public class BruteCheck {

    public static Random rnd;

    // ---------- random tree generation ----------
    static String randomTree(List<String> labels, double collapseProb) {
        List<String> parts = new ArrayList<>(labels);
        Collections.shuffle(parts, rnd);
        // random binary agglomeration, optionally merging 3+ at once (polytomy)
        while (parts.size() > 3) {
            int k = (rnd.nextDouble() < collapseProb && parts.size() > 4) ? 3 + rnd.nextInt(Math.min(3, parts.size() - 3)) : 2;
            List<String> chosen = new ArrayList<>();
            for (int i = 0; i < k; i++) chosen.add(parts.remove(rnd.nextInt(parts.size())));
            parts.add("(" + String.join(",", chosen) + ")");
        }
        return "(" + String.join(",", parts) + ");";
    }

    // ---------- leaf distances on the unrooted tree ----------
    static Map<Integer, int[]> leafDistances(Tree t, int nTaxa) {
        Map<TreeNode, List<TreeNode>> adj = new HashMap<>();
        for (TreeNode v : t.nodes) adj.put(v, new ArrayList<>());
        for (TreeNode v : t.nodes) {
            if (v.childs != null) for (TreeNode c : v.childs) { adj.get(v).add(c); adj.get(c).add(v); }
        }
        Map<Integer, int[]> res = new HashMap<>();
        for (TreeNode leaf : t.nodes) {
            if (!leaf.isLeaf()) continue;
            int[] d = new int[nTaxa];
            Arrays.fill(d, -1);
            Map<TreeNode, Integer> dist = new HashMap<>();
            ArrayDeque<TreeNode> q = new ArrayDeque<>();
            dist.put(leaf, 0); q.add(leaf);
            while (!q.isEmpty()) {
                TreeNode u = q.poll();
                if (u.isLeaf()) d[u.taxon.id] = dist.get(u);
                for (TreeNode w : adj.get(u)) if (!dist.containsKey(w)) { dist.put(w, dist.get(u) + 1); q.add(w); }
            }
            res.put(leaf.taxon.id, d);
        }
        return res;
    }

    // +1 if ab|cd, -1 if resolved otherwise, 0 if unresolved (four-point condition, unit lengths)
    static int topo(Map<Integer, int[]> D, int a, int b, int c, int d) {
        int s1 = D.get(a)[b] + D.get(c)[d];
        int s2 = D.get(a)[c] + D.get(b)[d];
        int s3 = D.get(a)[d] + D.get(b)[c];
        if (s1 < s2 && s1 < s3) return 1;
        if (s1 == s2 && s2 == s3) return 0;
        return -1;
    }

    static double bruteScore(List<Map<Integer, int[]>> dists, TaxaPerLevelWithPartition tx, int n) {
        double total = 0;
        for (Map<Integer, int[]> D : dists) {
            List<Integer> A = new ArrayList<>(), B = new ArrayList<>();
            for (int id : D.keySet()) {
                if (!tx.isInRealTaxa(id) && !tx.isInDummyTaxa(id)) continue;
                (tx.inWhichPartition(id) == 0 ? A : B).add(id);
            }
            for (int i = 0; i < A.size(); i++) for (int j = i + 1; j < A.size(); j++) {
                int a = A.get(i), b = A.get(j);
                if (sameDummy(tx, a, b)) continue;
                for (int k = 0; k < B.size(); k++) for (int l = k + 1; l < B.size(); l++) {
                    int c = B.get(k), d = B.get(l);
                    if (sameDummy(tx, c, d) || sameDummy(tx, a, c) || sameDummy(tx, a, d) || sameDummy(tx, b, c) || sameDummy(tx, b, d)) continue;
                    double w = tx.getWeight(a) * tx.getWeight(b) * tx.getWeight(c) * tx.getWeight(d);
                    total += w * topo(D, a, b, c, d);
                }
            }
        }
        return total;
    }

    static boolean sameDummy(TaxaPerLevelWithPartition tx, int a, int b) {
        return tx.isInDummyTaxa(a) && tx.isInDummyTaxa(b) && tx.inWhichDummyTaxa(a) == tx.inWhichDummyTaxa(b);
    }

    // ---------- random subproblem with nested dummies ----------
    static Object[] randomSubproblem(RealTaxon[] all) {
        List<RealTaxon> pool = new ArrayList<>(Arrays.asList(all));
        Collections.shuffle(pool, rnd);
        int nDummy = rnd.nextInt(Math.min(4, pool.size() / 3) + 1);
        List<DummyTaxon> dummies = new ArrayList<>();
        for (int i = 0; i < nDummy; i++) {
            int take = 1 + rnd.nextInt(3);
            List<RealTaxon> rs = new ArrayList<>();
            for (int j = 0; j < take && pool.size() > 4; j++) rs.add(pool.remove(pool.size() - 1));
            if (rs.isEmpty()) break;
            DummyTaxon[] nested = new DummyTaxon[0];
            if (!dummies.isEmpty() && rnd.nextBoolean()) {               // nest a previous dummy inside this one
                nested = new DummyTaxon[]{dummies.remove(rnd.nextInt(dummies.size()))};
            }
            dummies.add(new DummyTaxon(rs.toArray(new RealTaxon[0]), nested));
        }
        RealTaxon[] rts = pool.toArray(new RealTaxon[0]);
        DummyTaxon[] dts = dummies.toArray(new DummyTaxon[0]);
        int[] rp = new int[rts.length], dp = new int[dts.length];
        for (int i = 0; i < rp.length; i++) rp[i] = rnd.nextInt(2);
        for (int i = 0; i < dp.length; i++) dp[i] = rnd.nextInt(2);
        return new Object[]{rts, dts, rp, dp};
    }

    public static void main(String[] args) throws Exception {
        long seed = Long.parseLong(args[0]);
        int trials = Integer.parseInt(args[1]);
        String tmpDir = args[2];
        rnd = new Random(seed);
        if (args.length > 3 && args[3].equals("B")) src.Config.NON_QUARTET_TYPE = src.Config.NonQuartetType.B;
        double maxRelErrScore = 0, maxRelErrGain = 0;
        int checks = 0;

        for (int t = 0; t < trials; t++) {
            int n = 6 + rnd.nextInt(9);
            int k = 1 + rnd.nextInt(6);
            double collapse = rnd.nextDouble() < 0.5 ? 0 : 0.4;
            List<String> labels = new ArrayList<>();
            for (int i = 0; i < n; i++) labels.add("t" + i);
            File f = new File(tmpDir, "gt_" + seed + "_" + t + ".tre");
            try (PrintWriter pw = new PrintWriter(f)) {
                for (int g = 0; g < k; g++) {
                    List<String> present = new ArrayList<>(labels);
                    if (g > 0 && rnd.nextDouble() < 0.4) {               // incomplete gene tree (keep first tree complete so all taxa exist)
                        int drop = rnd.nextInt(Math.max(1, n - 5));
                        for (int i = 0; i < drop; i++) present.remove(rnd.nextInt(present.size()));
                    }
                    String tree = randomTree(present, collapse);
                    pw.println(tree);
                    if (rnd.nextDouble() < 0.3) pw.println(tree);           // duplicates exercise tripartition frequency dedup
                }
            }

            RealTaxon.count = 0;                                         // static id counter; must reset per dataset
            GeneTrees gts = new GeneTrees(f.getPath());
            gts.readTaxaNames();
            gts.readGeneTrees(null);
            List<Map<Integer, int[]>> dists = new ArrayList<>();
            for (Tree tr : gts.geneTrees) dists.add(leafDistances(tr, n));

            Object[] sp = randomSubproblem(gts.taxa);
            RealTaxon[] rts = (RealTaxon[]) sp[0];
            DummyTaxon[] dts = (DummyTaxon[]) sp[1];
            int[] rp = (int[]) sp[2], dp = (int[]) sp[3];
            if (rts.length + dts.length < 4) continue;
            TaxaPerLevelWithPartition tx = new TaxaPerLevelWithPartition(rts, dts, rp, dp, n);

            BookKeepingPerLevelv2 book = new BookKeepingPerLevelv2(gts);
            book.resetBookKeeping(tx);

            double fast = book.calculateScore();
            double brute = bruteScore(dists, tx, n);
            double rel = Math.abs(fast - brute) / Math.max(1, Math.abs(brute));
            maxRelErrScore = Math.max(maxRelErrScore, rel);
            checks++;
            if (rel > 1e-9) {
                System.out.printf("SCORE MISMATCH trial=%d fast=%.10f brute=%.10f file=%s%n", t, fast, brute, f);
            }

            // gains: compare against brute(after swap) - brute(before)
            double[][] rtG = new double[rts.length][2];
            double[] dtG = new double[dts.length];
            book.calculateScoreAndGains(rtG, dtG);
            for (int i = 0; i < rts.length; i++) {
                int p = tx.inWhichPartitionRealTaxonByIndex(i);
                tx.swapPartitionRealTaxon(i);
                double after = bruteScore(dists, tx, n);
                tx.swapPartitionRealTaxon(i);
                double expect = after - brute;
                double r = Math.abs(rtG[i][p] - expect) / Math.max(1, Math.abs(expect));
                maxRelErrGain = Math.max(maxRelErrGain, r);
                if (r > 1e-9) System.out.printf("RT GAIN MISMATCH trial=%d taxon=%s fast=%.8f brute=%.8f%n", t, rts[i], rtG[i][p], expect);
            }
            for (int i = 0; i < dts.length; i++) {
                tx.swapPartitionDummyTaxon(i);
                double after = bruteScore(dists, tx, n);
                tx.swapPartitionDummyTaxon(i);
                double expect = after - brute;
                double r = Math.abs(dtG[i] - expect) / Math.max(1, Math.abs(expect));
                maxRelErrGain = Math.max(maxRelErrGain, r);
                if (r > 1e-9) System.out.printf("DT GAIN MISMATCH trial=%d dummy=%d fast=%.8f brute=%.8f%n", t, i, dtG[i], expect);
            }

            // incremental swap path: actually swap a taxon through the book, then recheck score
            if (rts.length > 0) {
                int i = rnd.nextInt(rts.length);
                book.swapTaxon(i, false);
                double fast2 = book.calculateScore();
                double brute2 = bruteScore(dists, tx, n);
                double r2 = Math.abs(fast2 - brute2) / Math.max(1, Math.abs(brute2));
                maxRelErrScore = Math.max(maxRelErrScore, r2);
                if (r2 > 1e-9) System.out.printf("POST-RT-SWAP MISMATCH trial=%d fast=%.8f brute=%.8f%n", t, fast2, brute2);
            }
            if (dts.length > 0) {
                int i = rnd.nextInt(dts.length);
                book.swapTaxon(i, true);
                double fast3 = book.calculateScore();
                double brute3 = bruteScore(dists, tx, n);
                double r3 = Math.abs(fast3 - brute3) / Math.max(1, Math.abs(brute3));
                maxRelErrScore = Math.max(maxRelErrScore, r3);
                if (r3 > 1e-9) System.out.printf("POST-DT-SWAP MISMATCH trial=%d fast=%.8f brute=%.8f%n", t, fast3, brute3);
            }
            f.delete();
        }
        System.out.printf("seed=%d subproblems=%d  max rel err score=%.2e  gains=%.2e%n", seed, checks, maxRelErrScore, maxRelErrGain);
        if (maxRelErrScore > 1e-9 || maxRelErrGain > 1e-9) System.exit(1);
    }
}
