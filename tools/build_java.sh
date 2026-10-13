#!/usr/bin/env bash
# Build the Java side used for validation and measurement, from a COPY of the wQFM-TREE sources.
#
#   build/java-ref    original algorithm + TopLevelBench harness + brute-force checkers
#   build/java-instr  same algorithm with Phase-0 measurement output (tools/java/instrument.py)
#
# The upstream tree does not compile as-is: the legacy v1 classes no longer match the current
# interfaces, and QFM2 has an unused import that drags them in. Both are excluded here.
#
# usage: tools/build_java.sh [wQFM-TREE dir, default ..]
set -euo pipefail

HERE="$(cd "$(dirname "$0")/.." && pwd)"
WQ="$(cd "${1:-$HERE/..}" && pwd)"
B="$HERE/build"
LEGACY='src/(QFM|DebugMain|DSPerLevel/BookKeepingPerLevel|DSPerLevel/BookKeepingPerTree|ScoreCalculator/NumSatCalculatorBinaryNode|ScoreCalculator/NumSatCalculatorNodeE)\.java$'

[ -d "$WQ/src" ] || { echo "no src/ in $WQ" >&2; exit 1; }

prepare() {
    local dir=$1
    rm -rf "$dir"
    mkdir -p "$dir"
    cp -R "$WQ/src" "$dir/src"
    find "$dir" -name '._*' -delete
    sed -i.bak '/^import src\.DSPerLevel\.BookKeepingPerLevel;$/d' "$dir/src/QFM2.java" && rm "$dir/src/QFM2.java.bak"
    cp "$HERE/tools/java/Phase0Stats.java" "$dir/src/"
}

compile() {
    local dir=$1
    # ._* AppleDouble files reappear on exFAT whenever a file is written, so filter them here, not just delete
    (cd "$dir" && find src -name '*.java' ! -name '._*' | grep -v -E "$LEGACY" > sources.txt \
        && mkdir -p classes && javac -nowarn -sourcepath "" -d classes @sources.txt)
}

prepare "$B/java-ref"
cp "$HERE/tools/java/TopLevelBench.java" "$B/java-ref/src/"
compile "$B/java-ref"
javac -nowarn -cp "$B/java-ref/classes" -d "$B/java-ref/classes" \
    "$HERE/tools/java/BruteCheck.java" "$HERE/tools/java/CompressCheck.java"

prepare "$B/java-instr"
python3 "$HERE/tools/java/instrument.py" "$B/java-instr"
compile "$B/java-instr"

echo "built $B/java-ref/classes and $B/java-instr/classes from $WQ/src"
