#!/usr/bin/env bash
# Phase-0 measurement on one gene-tree file:
#   clean trees -> greedy consensus (or a given one) -> instrumented wQFM-TREE -> summary.
#
# usage: tools/phase0/run_phase0.sh <gene trees> <out dir> [consensus tree]
#   env JAVA_XMX (default 8g)
# Pass a consensus tree for incomplete gene trees (the built-in greedy consensus needs complete trees),
# or PAUP's greedy consensus to reproduce run.sh exactly.
set -euo pipefail

HERE="$(cd "$(dirname "$0")/../.." && pwd)"
gt=$1
out=$2
cons=${3:-}
[ -d "$HERE/build/java-instr/classes" ] || { echo "run 'make java' first" >&2; exit 1; }
mkdir -p "$out"

python3 -I "$HERE/tools/greedy_consensus.py" --clean "$gt" > "$out/gt.tre"
if [ -n "$cons" ]; then
    python3 -I "$HERE/tools/greedy_consensus.py" --clean "$cons" > "$out/cons.tre"
else
    python3 -I "$HERE/tools/greedy_consensus.py" "$out/gt.tre" > "$out/cons.tre"
fi

start=$(date +%s)
java -Xmx"${JAVA_XMX:-8g}" -cp "$HERE/build/java-instr/classes" src.Main \
    "$out/gt.tre" "$out/cons.tre" "$out/species.tre" A > "$out/run.log" 2> "$out/run.err"
echo "instrumented wQFM-TREE finished in $(( $(date +%s) - start )) s" | tee "$out/summary.txt"
python3 -I "$HERE/tools/phase0/analyze.py" "$out/run.log" | tee -a "$out/summary.txt"
