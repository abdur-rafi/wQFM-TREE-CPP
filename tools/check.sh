#!/usr/bin/env bash
# Correctness checks; each step exits non-zero on any mismatch.
#  1. BruteCheck:    Java score and gains (real + dummy taxa, after swaps) vs explicit quartet enumeration
#  2. CompressCheck: collapsing single-dummy subtrees preserves the score exactly; Phase0Stats counts are right
#  3. compare.sh:    C++ prototype vs Java on random incomplete gene trees with polytomies (bit-identical)
#  4. brute_check:   C++ scorer (score, real and dummy gains, swaps) vs explicit quartet enumeration
#  5. validate.sh:   full search, C++ replay of a Java trace: same tree, every score and gain within 1e-9
#  6. threads:       byte-identical output with 1 and 4 threads
#  7. consensus:     built-in greedy consensus equals PAUP's on the repo's 11-taxon example
set -euo pipefail

HERE="$(cd "$(dirname "$0")/.." && pwd)"
B="$HERE/build"
mkdir -p "$B/tmp"
quiet() { grep -v -E '^(taxon count|Gene trees count|total internal nodes|unique partitions)' || true; }

echo "== C++ scorer vs brute force (incl. subproblems with no real taxa)"
"$B/brute_check" 1 300
"$B/brute_check" 2 300
echo "== BruteCheck (Java scorer)"
java -cp "$B/java-ref/classes" BruteCheck 1 300 "$B/tmp" | quiet
echo "== CompressCheck"
java -cp "$B/java-ref/classes" CompressCheck 11 300 "$B/tmp" | quiet

echo "== C++ vs Java (incomplete trees with polytomies)"
python3 - "$B/tmp/inc_poly.tre" <<'EOF'
import random, sys
rnd = random.Random(7)
n, k = 60, 40
def rtree(labels):
    parts = list(labels)
    rnd.shuffle(parts)
    while len(parts) > 3:
        m = 3 if rnd.random() < 0.25 and len(parts) > 4 else 2
        kids = [parts.pop(rnd.randrange(len(parts))) for _ in range(m)]
        parts.append('(' + ','.join(kids) + ')')
    return '(' + ','.join(parts) + ');'
with open(sys.argv[1], 'w') as f:
    labels = [f't{i}' for i in range(n)]
    f.write(rtree(labels) + '\n')                      # one complete tree so every taxon appears
    for _ in range(k - 1):
        f.write(rtree([t for t in labels if rnd.random() > 0.15]) + '\n')
EOF
"$HERE/tools/compare/compare.sh" "$B/tmp/inc_poly.tre" "$B/tmp/compare" 4 | grep -E 'identical|RESULT'
echo "== full search: Java trace vs C++ replay (incomplete gene trees with polytomies)"
(cd "$HERE" && python3 tools/sim/gen_only.py 40 40 6 5 "$B/tmp/search.tre" --missing 0.15 --collapse 0.25 --species "$B/tmp/search.sp.tre" > /dev/null)
"$HERE/tools/compare/validate.sh" "$B/tmp/search.tre" "$B/tmp/search.sp.tre" "$B/tmp/validate" | grep -E 'RF|RESULT'
echo "== threads: identical tree with 1 and 4 threads"
"$B/wqfm-tree" -i "$B/tmp/search.tre" -c "$B/tmp/search.sp.tre" -o "$B/tmp/t1.tre" -t 1 2> /dev/null
"$B/wqfm-tree" -i "$B/tmp/search.tre" -c "$B/tmp/search.sp.tre" -o "$B/tmp/t4.tre" -t 4 2> /dev/null
cmp "$B/tmp/t1.tre" "$B/tmp/t4.tre" && echo "identical"
echo "== built-in greedy consensus vs PAUP (11-taxon example in ../input)"
"$B/wqfm-tree" -i "$HERE/../input/gtree_11tax_est_5genes_R1.tre-cleaned" -o "$B/tmp/c11.tre" --write-consensus "$B/tmp/c11.cons" 2> /dev/null
python3 -I "$HERE/tools/compare/rf.py" "$HERE/../input/gtree_11tax_est_5genes_R1.tre-cleaned-cons.greedy.tree" "$B/tmp/c11.cons"
echo "all checks passed"
