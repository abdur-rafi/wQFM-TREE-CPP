#!/usr/bin/env bash
# Validate the C++ version against Java on ASTRAL-II replicates: both get the same gene trees (the
# first k of the replicate) and the same consensus (the built-in one), then tools/compare/validate.sh
# runs traced Java, an independent C++ run and a C++ replay of Java's decisions.
#
# usage: tools/astral2/validate_batch.sh <astral2 data dir> <out dir> <condition:rep:genes>...
#   e.g. tools/astral2/validate_batch.sh data/astral2 runs/validate model.200.500000.0.000001:01:50
set -uo pipefail

HERE="$(cd "$(dirname "$0")/../.." && pwd)"
data=$1
out=$2
shift 2
fail=0
for spec in "$@"; do
    IFS=: read -r cond rep genes <<< "$spec"
    o="$out/$cond-$rep-$genes"
    mkdir -p "$o"
    echo "=== $cond, replicate $rep, $genes genes"
    grep '(' "$data/$cond/$rep/estimatedgenetre" | head -n "$genes" > "$o/gt.tre"
    "$HERE/build/wqfm-tree" -i "$o/gt.tre" -o "$o/native.tre" --write-consensus "$o/cons.tre" 2> "$o/native.log" \
        || { echo "C++ run failed, see $o/native.log"; fail=1; continue; }
    SECONDS=0
    "$HERE/tools/compare/validate.sh" "$o/gt.tre" "$o/cons.tre" "$o/v"
    status=$?
    echo "validate.sh exit $status after $SECONDS s"
    [ $status -eq 0 ] || fail=1
done
exit $fail
