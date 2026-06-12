#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

make

mkdir -p results figures

./counting \
  1.0 \
  graphs/facebook/facebook.txt \
  graphs/facebook/facebook_partitioned_80 \
  results/facebook_metrics.csv | tee results/facebook_run.log

python3 scripts/plot_facebook_metrics.py \
  --input results/facebook_metrics.csv \
  --png figures/facebook_metrics.png \
  --pdf figures/facebook_metrics.pdf

echo
echo "Saved:"
echo "  results/facebook_run.log"
echo "  results/facebook_metrics.csv"
echo "  figures/facebook_metrics.png"
echo "  figures/facebook_metrics.pdf"
