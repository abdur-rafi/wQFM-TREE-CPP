"""Add step-by-step trace output (tools/java/Trace.java) to a COPY of the wQFM-TREE Java source.

usage: trace_instrument.py <dir containing src/>

Enabled at run time with -Dwqfm.trace=<file>; without it the build behaves exactly like the original.
Each anchor must occur exactly once, so a changed upstream source fails loudly.
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


patch(root + '/src/Main.java', [
    ("        trees.readGeneTrees(consensusTreePartition.dist);\n",
     "        trees.readGeneTrees(consensusTreePartition.dist);\n        Trace.taxa(trees.taxa);\n"),
])

patch(root + '/src/QFM2.java', [
    ("        this.level++;\n",
     "        this.level++;\n        Trace.sub(book.taxas);\n"),
    ("        var x = book.divide(initPartition, allowSingleton);\n",
     "        Trace.fin(book.taxas);\n        var x = book.divide(initPartition, allowSingleton);\n"),
    ("                cg += gain;\n",
     "                cg += gain;\n                Trace.move(x.isDummy, x.index, gain, cg);\n"),
    ("        if(maxCgIndex == -1){\n",
     "        Trace.end(maxCgIndex);\n        if(maxCgIndex == -1){\n"),
])

patch(root + '/src/InitialPartition/ConsensusTreePartition.java', [
    ("    double scoreForPartitionByNode(TreeNode node, RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton){\n",
     "    double scoreForPartitionByNode(TreeNode node, RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton){\n"
     "        return src.Trace.cand(scoreForPartitionByNodeInner(node, rts, dts, allowSingleton));\n"
     "    }\n\n"
     "    double scoreForPartitionByNodeInner(TreeNode node, RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton){\n"),
    ("    @Override\n    public MakePartitionReturnType makePartition(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {\n",
     "    @Override\n    public MakePartitionReturnType makePartition(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {\n"
     "        src.Trace.part(rts, dts);\n"
     "        MakePartitionReturnType res = makePartitionInner(rts, dts, allowSingleton);\n"
     "        src.Trace.result(res);\n"
     "        return res;\n"
     "    }\n\n"
     "    public MakePartitionReturnType makePartitionInner(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {\n"),
])
print('trace instrumentation applied')
