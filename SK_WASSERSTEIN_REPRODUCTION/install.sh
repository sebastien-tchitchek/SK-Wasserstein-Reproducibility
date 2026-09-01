#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE1="$ROOT/01_COMPUTE_DIAGRAMS_AND_MATRICES"
STAGE2="$ROOT/02_REPRODUCE_FIGURES_AND_TABLES"
RESULTS="$ROOT/results"
LOG="$ROOT/.run.log"
START="$(date +%s)"
mkdir -p "$RESULTS"
: > "$LOG"
format_time() {
  local seconds="$1"
  printf '%02d:%02d:%02d' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
}
fail_run() {
  local code="$1"
  local elapsed="$(( $(date +%s) - START ))"
  printf 'Computation failed after %s.\n' "$(format_time "$elapsed")" >&2
  printf 'Log: %s\n' "$LOG" >&2
  exit "$code"
}
if bash "$STAGE1/run.sh" >> "$LOG" 2>&1; then
  :
else
  code="$?"
  fail_run "$code"
fi
if bash "$STAGE2/run.sh" >> "$LOG" 2>&1; then
  :
else
  code="$?"
  fail_run "$code"
fi
rm -rf "$RESULTS/normalized_diagrams" "$RESULTS/matrices"
ln -s "../01_COMPUTE_DIAGRAMS_AND_MATRICES/NORMALIZED_DIAGRAM_RESULTS" "$RESULTS/normalized_diagrams"
ln -s "../01_COMPUTE_DIAGRAMS_AND_MATRICES/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS" "$RESULTS/matrices"
rm -f "$STAGE1/01_diagrams.log" "$STAGE1/02_matrices.log" "$STAGE1/ttk_configuration.sh"
rm -f "$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS.zip"
rm -rf "$STAGE1/.cache_wmt_reproduction"
rm -f "$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS/README.txt" "$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS/SHA256SUMS.txt"
rm -f "$STAGE2/analysis.log"
rm -rf "$STAGE2/.venv_reproduction" "$STAGE2/.extracted_matrices"
rm -f "$LOG"
ELAPSED="$(( $(date +%s) - START ))"
printf '%s\n' "$(format_time "$ELAPSED")" > "$RESULTS/computation_time.txt"
printf 'Computation time: %s\n' "$(format_time "$ELAPSED")"
printf 'Results: %s\n' "$RESULTS"
