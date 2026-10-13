package src;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;

import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.InitialPartition.MakePartitionReturnType;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;

/**
 * Step-by-step trace of the search, enabled with -Dwqfm.trace=<file>. The format is shared with the
 * C++ version (src/trace.hpp there); see that file for the line types. Read-only: the search is unchanged.
 */
public final class Trace {
    private static final PrintWriter out = open();

    private static PrintWriter open() {
        String path = System.getProperty("wqfm.trace");
        if (path == null) return null;
        try {
            PrintWriter w = new PrintWriter(new BufferedWriter(new FileWriter(path), 1 << 16));
            Runtime.getRuntime().addShutdownHook(new Thread(w::close));
            return w;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void taxa(RealTaxon[] taxa) {
        if (out == null) return;
        for (RealTaxon t : taxa) out.println("T " + t.id + " " + t.label);
    }

    static String desc(DummyTaxon d) {
        StringBuilder sb = new StringBuilder("{");
        for (int i = 0; i < d.realTaxa.length; i++) sb.append(i > 0 ? "," : "").append(d.realTaxa[i].id);
        sb.append('|');
        for (int i = 0; i < d.dummyTaxa.length; i++) sb.append(i > 0 ? ";" : "").append(desc(d.dummyTaxa[i]));
        return sb.append('}').toString();
    }

    static String key(RealTaxon[] rts, DummyTaxon[] dts) {
        StringBuilder sb = new StringBuilder("R:");
        for (int i = 0; i < rts.length; i++) sb.append(i > 0 ? "," : "").append(rts[i].id);
        sb.append(" D:");
        for (int i = 0; i < dts.length; i++) sb.append(i > 0 ? ";" : "").append(desc(dts[i]));
        return sb.toString();
    }

    static String bits(int[] r, int[] d) {
        StringBuilder sb = new StringBuilder();
        for (int x : r) sb.append(x);
        sb.append('|');
        for (int x : d) sb.append(x);
        return sb.toString();
    }

    public static void part(RealTaxon[] rts, DummyTaxon[] dts) { if (out != null) out.println("P " + key(rts, dts)); }

    public static double cand(double score) { if (out != null) out.println("C " + score); return score; }

    public static void result(MakePartitionReturnType r) {
        if (out != null) out.println("R " + bits(r.realTaxonPartition, r.dummyTaxonPartition));
    }

    public static void sub(TaxaPerLevelWithPartition t) {
        if (out == null) return;
        out.println("S " + key(t.realTaxa, t.dummyTaxa));
        out.println("I " + bits(t.realTaxonPartition, t.dummyTaxonPartition));
    }

    public static void move(boolean dummy, int index, double gain, double cg) {
        if (out != null) out.println("M " + (dummy ? 1 : 0) + " " + index + " " + gain + " " + cg);
    }

    public static void end(int maxCgIndex) { if (out != null) out.println("E " + maxCgIndex); }

    public static void fin(TaxaPerLevelWithPartition t) {
        if (out != null) out.println("F " + bits(t.realTaxonPartition, t.dummyTaxonPartition));
    }
}
