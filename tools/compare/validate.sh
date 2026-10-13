#!/usr/bin/env bash
# Full-search validation of the C++ version against the Java version on one dataset:
#   1. traced Java run (build/java-trace)
#   2. independent C++ run with Java's taxon order       -> agreement until the first tie broken differently
#   3. C++ replay of Java's decisions                     -> every score and gain compared; same tree expected
#
# usage: tools/compare/validate.sh <gene trees> <consensus tree> <out dir>
#   env JAVA_XMX (default 6g)
set -uo pipefail

HERE="$(cd "$(dirname "$0")/../.." && pwd)"
gt=$1
cons=$2
out=$3
[ -x "$HERE/build/wqfm-tree" ] && [ -d "$HERE/build/java-trace/classes" ] || { echo "run 'make && make java' first" >&2; exit 1; }
mkdir -p "$out"
python3 -I "$HERE/tools/greedy_consensus.py" --clean "$gt" > "$out/gt.tre"
python3 -I "$HERE/tools/greedy_consensus.py" --clean "$cons" > "$out/cons.tre"

secs() { python3 -c 'import time; print(time.time())'; }

t0=$(secs)
java -Xmx"${JAVA_XMX:-6g}" -Dwqfm.trace="$out/java.trace" -cp "$HERE/build/java-trace/classes" src.Main \
    "$out/gt.tre" "$out/cons.tre" "$out/java.tre" A > "$out/java.log" 2>&1 || { echo "Java run failed, see $out/java.log"; exit 1; }
t1=$(secs)
"$HERE/build/wqfm-tree" -i "$out/gt.tre" -c "$out/cons.tre" -o "$out/cpp.tre" \
    --taxa-order "$out/java.trace" --trace "$out/cpp.trace" 2> "$out/cpp.log" || { cat "$out/cpp.log"; exit 1; }
t2=$(secs)
"$HERE/build/wqfm-tree" -i "$out/gt.tre" -c "$out/cons.tre" -o "$out/replay.tre" \
    --replay "$out/java.trace" --trace "$out/replay.trace" 2> "$out/replay.log" || { cat "$out/replay.log"; exit 1; }

python3 - "$t0" "$t1" "$t2" <<'EOF'
import sys
t0, t1, t2 = map(float, sys.argv[1:])
print(f"time (traced runs): Java {t1 - t0:.1f} s, C++ {t2 - t1:.1f} s")
EOF
echo "== independent C++ run vs Java"
python3 -I "$HERE/tools/compare/rf.py" "$out/java.tre" "$out/cpp.tre"
python3 -I "$HERE/tools/compare/trace_compare.py" "$out/java.trace" "$out/cpp.trace" | sed 's/^/   /'
echo "== C++ replay of Java's decisions"
python3 -I "$HERE/tools/compare/rf.py" "$out/java.tre" "$out/replay.tre"
rf_status=$?
python3 -I "$HERE/tools/compare/trace_compare.py" "$out/java.trace" "$out/replay.trace" | sed 's/^/   /'
cmp_status=$?
# replay must reproduce Java's tree and agree on every number: that is the pass/fail criterion
exit $(( rf_status != 0 || cmp_status != 0 ))
