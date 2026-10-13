#!/usr/bin/env bash
# Correctness checks; each step exits non-zero on any mismatch.
#  1. BruteCheck:    Java score and gains (real + dummy taxa, after swaps) vs explicit quartet enumeration
#  2. CompressCheck: collapsing single-dummy subtrees preserves the score exactly; Phase0Stats counts are right
#  3. compare.sh:    C++ prototype vs Java on random incomplete gene trees with polytomies (bit-identical)
set -euo pipefail

HERE="$(cd "$(dirname "$0")/.." && pwd)"
B="$HERE/build"
mkdir -p "$B/tmp"
quiet() { grep -v -E '^(taxon count|Gene trees count|total internal nodes|unique partitions)' || true; }

echo "== BruteCheck"
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
echo "all checks passed"
