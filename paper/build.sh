#!/bin/bash
# Build main.pdf with tectonic. Its package cache goes to the external SSD (the internal disk is
# nearly full); override with TECTONIC_CACHE_DIR.
set -euo pipefail
cd "$(dirname "$0")"
export TECTONIC_CACHE_DIR="${TECTONIC_CACHE_DIR:-/Volumes/Extreme SSD/research/tectonic-cache}"
mkdir -p "$TECTONIC_CACHE_DIR"
tectonic --keep-logs main.tex
