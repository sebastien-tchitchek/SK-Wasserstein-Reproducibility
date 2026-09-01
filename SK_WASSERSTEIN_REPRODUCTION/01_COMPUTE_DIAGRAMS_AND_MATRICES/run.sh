#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================"
echo "STAGE 1/2 — NORMALIZED PERSISTENCE DIAGRAMS"
echo "============================================================"
bash "$ROOT/scripts/diagrams/run_diagrams.sh"

echo
echo "============================================================"
echo "STAGE 2/2 — SK WASSERSTEIN, W_GAMMA, AND W2 MATRICES"
echo "============================================================"
bash "$ROOT/scripts/matrices/run_matrices.sh"
