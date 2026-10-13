#!/usr/bin/env bash
# Java vs C++ on the same top-level partition: one gain sweep and one full FM pass.
# Reports whether scores, gains, every FM move and the final partition agree exactly, plus timings.
#
# usage: tools/compare/compare.sh <gene trees> <out dir> [threads, default 4]
#   env JAVA_XMX (default 6g), SEED (partition seed, default 42)
set -euo pipefail

HERE="$(cd "$(dirname "$0")/../.." && pwd)"
gt=$1
out=$2
T=${3:-4}
[ -x "$HERE/build/gainsweep" ] && [ -d "$HERE/build/java-ref/classes" ] || { echo "run 'make && make java' first" >&2; exit 1; }
mkdir -p "$out"
p="$out/run"
python3 -I "$HERE/tools/greedy_consensus.py" --clean "$gt" > "$out/gt.tre"

java_bench() { java -Xmx"${JAVA_XMX:-6g}" -cp "$HERE/build/java-ref/classes" src.TopLevelBench "$out/gt.tre" "$p" "${SEED:-42}" "$@"; }

echo "== gain sweep (Java, then C++ with $T threads)"
java_bench 10 sweep | grep -E '^(score|sweep_ms)'
"$HERE/build/gainsweep" "$out/gt.tre" "$p.taxa" "$p.part" "$T" 50 sweep "$p" | grep -E '^(score|sweep_ms|ns_per_node)'
echo "== one FM pass (Java, then C++)"
java_bench 0 fm | grep -E '^fm_'
"$HERE/build/gainsweep" "$out/gt.tre" "$p.taxa" "$p.part" "$T" 0 fm "$p" | grep -E '^fm_'
echo "== agreement"
python3 -I "$HERE/tools/compare/compare.py" "$p"
