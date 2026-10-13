package src;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

import src.DSPerLevel.BookKeepingPerLevelv2;
import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.PreProcessing.GeneTrees;
import src.Taxon.DummyTaxon;

/** Top-level (no dummy) gain-sweep benchmark and FM-pass trace for comparison with the C++ prototype. */
public class TopLevelBench {

    public static void main(String[] args) throws Exception {
        String gtPath = args[0], out = args[1], mode = args[4];
        long seed = Long.parseLong(args[2]);
        int reps = Integer.parseInt(args[3]);

        GeneTrees trees = new GeneTrees(gtPath);
        trees.readTaxaNames();
        trees.readGeneTrees(null);
        int n = trees.taxa.length;

        int[] part = new int[n];
        Random rnd = new Random(seed);
        int c1;
        do { c1 = 0; for (int i = 0; i < n; i++) { part[i] = rnd.nextInt(2); c1 += part[i]; } } while (c1 < 3 || n - c1 < 3);
        try (PrintWriter pw = new PrintWriter(out + ".taxa")) { for (int i = 0; i < n; i++) pw.println(i + " " + trees.taxa[i].label); }
        try (PrintWriter pw = new PrintWriter(out + ".part")) { for (int i = 0; i < n; i++) pw.println(i + " " + part[i]); }

        TaxaPerLevelWithPartition x = new TaxaPerLevelWithPartition(trees.taxa, new DummyTaxon[0], part.clone(), new int[0], n);
        BookKeepingPerLevelv2 book = new BookKeepingPerLevelv2(trees);
        book.resetBookKeeping(x);

        if (mode.equals("sweep")) {
            double[][] g = null;
            double sc = 0;
            for (int w = 0; w < 3; w++) { g = new double[n][2]; sc = book.calculateScoreAndGains(g, new double[0]); }
            long t0 = System.nanoTime();
            for (int r = 0; r < reps; r++) { g = new double[n][2]; sc = book.calculateScoreAndGains(g, new double[0]); }
            double ms = (System.nanoTime() - t0) / 1e6 / reps;
            System.out.printf("score %.17g%nsweep_ms %.4f%n", sc, ms);
            try (PrintWriter pw = new PrintWriter(out + ".java.gains")) {
                for (int i = 0; i < n; i++) pw.printf("%d %.17g%n", i, g[i][x.inWhichPartitionRealTaxonByIndex(i)]);
            }
            return;
        }

        // FM pass: QFM2.oneInteration with a move log
        long t0 = System.nanoTime();
        double EPS = 1e-5, cg = 0, maxCg = 0;
        int maxCgIndex = -1;
        boolean singleton = book.taxas.getTaxonCountInPartition(0) == 1 || book.taxas.getTaxonCountInPartition(1) == 1;
        boolean[] rtLocked = new boolean[n];
        boolean[] dtLocked = new boolean[0];
        ArrayList<QFM2.Swap> swaps = new ArrayList<>();
        try (PrintWriter log = new PrintWriter(out + ".java.fmlog")) {
            while (true) {
                double[][] rtGains = new double[n][2];
                book.calculateScoreAndGains(rtGains, new double[0]);
                QFM2.Swap s = QFM2.swapMax(book, rtGains, new double[0], rtLocked, dtLocked);
                if (s == null) break;
                swaps.add(s);
                cg += s.gain;
                log.printf("%d %.17g %.17g%n", s.index, s.gain, cg);
                if (singleton && maxCgIndex == -1) { maxCg = cg; maxCgIndex = swaps.size() - 1; }
                if (cg > maxCg && Math.abs(maxCg - cg) > EPS) { maxCg = cg; maxCgIndex = swaps.size() - 1; }
            }
        }
        if (maxCgIndex == -1) {
            if (swaps.size() != n) for (int i = swaps.size() - 1; i >= 0; --i) book.swapTaxon(swaps.get(i).index, false);
        } else {
            for (int i = swaps.size() - 1; i > maxCgIndex; --i) book.swapTaxon(swaps.get(i).index, false);
        }
        double secs = (System.nanoTime() - t0) / 1e9;
        double sc = book.calculateScore();
        System.out.printf("fm_moves %d kept %d improved %d max_cg %.17g final_score %.17g%nfm_pass_s %.3f%n",
            swaps.size(), maxCgIndex + 1, maxCgIndex == -1 ? 0 : 1, maxCg, sc, secs);
        try (PrintWriter pw = new PrintWriter(out + ".java.fmpart")) {
            for (int i = 0; i < n; i++) pw.println(i + " " + book.taxas.inWhichPartitionRealTaxonByIndex(i));
        }
    }
}
