"""Add Phase-0 measurement output to a COPY of the wQFM-TREE Java source (never the original).

usage: instrument.py <dir containing src/>

Read-only instrumentation: prints one P0SUB line per subproblem (compression statistics from
Phase0Stats, FM time, passes) and one P0CAND line per initial-split computation (candidates
considered / rescored / distinct, time). The inferred species tree is unchanged.
Each anchor string must occur exactly once, so a changed upstream source fails loudly.
"""
import sys
root = sys.argv[1]
def patch(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        n = s.count(old)
        assert n == 1, (path, n, old[:60])
        s = s.replace(old, new)
    open(path, 'w').write(s)

patch(root + '/src/QFM2.java', [
 ("        this.level++;\n        int itrCount = 0;\n",
  "        this.level++;\n        int itrCount = 0;\n"
  "        String p0stats = Phase0Stats.subproblem(geneTrees, book.taxas);\n"
  "        int p0s = book.taxas.realTaxonCount + book.taxas.dummyTaxonCount, p0r = book.taxas.realTaxonCount, p0d = book.taxas.dummyTaxonCount;\n"
  "        long p0t0 = System.nanoTime();\n"),
 ("        var x = book.divide(initPartition, allowSingleton);\n",
  "        long p0t1 = System.nanoTime();\n"
  "        var x = book.divide(initPartition, allowSingleton);\n"
  "        long p0t2 = System.nanoTime();\n"
  "        System.out.println(\"P0SUB level=\" + level + \" s=\" + p0s + \" r=\" + p0r + \" d=\" + p0d + \" passes=\" + (itrCount + 1)"
  " + \" fmMs=\" + (p0t1 - p0t0) / 1000000 + \" divideMs=\" + (p0t2 - p0t1) / 1000000 + \" \" + p0stats);\n"),
])

patch(root + '/src/InitialPartition/ConsensusTreePartition.java', [
 ("    double score;\n",
  "    double score;\n"
  "    long p0considered, p0rescored;\n"
  "    java.util.HashSet<Long> p0distinct = new java.util.HashSet<>();\n"
  "    static long p0mix(long z) { z += 0x9E3779B97F4A7C15L; z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L; z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL; return z ^ (z >>> 31); }\n"),
 ("        for(i = 0; i < dts.length; ++i){\n            if(node.info.branches[0].dummyTaxaWeightsIndividual[i] >= .5){\n                dtsp[i] = 1;\n            }\n        }\n\n        if(this.book == null){",
  "        for(i = 0; i < dts.length; ++i){\n            if(node.info.branches[0].dummyTaxaWeightsIndividual[i] >= .5){\n                dtsp[i] = 1;\n            }\n        }\n"
  "        p0considered++;\n"
  "        { long h = 0, hAll = 0;\n"
  "          for (int q = 0; q < rts.length; ++q) { long r = p0mix(2L * rts[q].id); hAll ^= r; if (rtsP[q] == 1) h ^= r; }\n"
  "          for (int q = 0; q < dts.length; ++q) { long r = p0mix(2L * dts[q].id + 1); hAll ^= r; if (dtsp[q] == 1) h ^= r; }\n"
  "          p0distinct.add(Math.min(h, h ^ hAll)); }\n"
  "\n        if(this.book == null){"),
 ("        this.score = this.book.calculateScore();\n",
  "        p0rescored++;\n        this.score = this.book.calculateScore();\n"),
 ("    @Override\n    public MakePartitionReturnType makePartition(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {\n",
  "    @Override\n    public MakePartitionReturnType makePartition(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {\n"
  "        p0considered = 0; p0rescored = 0; p0distinct.clear();\n"
  "        long p0t0 = System.nanoTime();\n"
  "        var p0res = makePartitionInner(rts, dts, allowSingleton);\n"
  "        System.out.println(\"P0CAND s=\" + (rts.length + dts.length) + \" r=\" + rts.length + \" d=\" + dts.length + \" considered=\" + p0considered"
  " + \" rescored=\" + p0rescored + \" distinct=\" + p0distinct.size() + \" ms=\" + (System.nanoTime() - p0t0) / 1000000);\n"
  "        return p0res;\n    }\n\n"
  "    public MakePartitionReturnType makePartitionInner(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {\n"),
])
print("patched OK")
