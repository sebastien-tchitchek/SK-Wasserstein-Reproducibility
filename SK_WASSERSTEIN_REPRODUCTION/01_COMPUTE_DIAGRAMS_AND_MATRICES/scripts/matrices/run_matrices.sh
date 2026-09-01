#!/usr/bin/env bash
set -Eeuo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$MODULE_DIR/../.." && pwd)"
DIAGRAMS_ROOT="$PROJECT_ROOT/NORMALIZED_DIAGRAM_RESULTS"
OUTPUT_ROOT="$PROJECT_ROOT/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS"
ARCHIVE="$PROJECT_ROOT/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS.zip"
CONFIG_FILE="$PROJECT_ROOT/ttk_configuration.sh"
LOG_FILE="$PROJECT_ROOT/02_matrices.log"
TIMING_RUNS="${TIMING_RUNS:-1}"
HEARTBEAT_SECONDS="${HEARTBEAT_SECONDS:-120}"
COLLECTION_TIMEOUT_MINUTES="${COLLECTION_TIMEOUT_MINUTES:-600}"
DELTA_LIM="${DELTA_LIM:-0.01}"

mkdir -p "$PROJECT_ROOT"
touch "$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

fail() { echo "ERROR: $*" >&2; echo "Log: $LOG_FILE" >&2; exit 1; }
on_error() {
  local code=$?
  echo
  echo "Matrix computation stopped. Rerun: bash run.sh"
  echo "Completed matrices and timing measurements will be reused."
  echo "Log: $LOG_FILE"
  exit "$code"
}
trap on_error ERR

[ -d "$DIAGRAMS_ROOT" ] || fail "Persistence diagrams are missing: $DIAGRAMS_ROOT"
[ -f "$DIAGRAMS_ROOT/2004_isabel_3D/diagrams.vtm" ] || fail "The persistence-diagram benchmark is incomplete."
[[ "$TIMING_RUNS" =~ ^[1-9][0-9]*$ ]] || fail "TIMING_RUNS must be an integer >= 1."

PVPYTHON_PATH="${PVPYTHON_PATH:-}"
TTK_PLUGIN_PATH="${TTK_PLUGIN_PATH:-}"
if [ -f "$CONFIG_FILE" ]; then
  source "$CONFIG_FILE"
fi

verify_candidate() {
  local pv="$1" plugin="${2:-}"
  [ -x "$pv" ] || return 1
  if [ -n "$plugin" ]; then
    [ -f "$plugin" ] || return 1
    TTK_PLUGIN="$plugin" "$pv" "$MODULE_DIR/verify_ttk_matrices.py" \
      --input "$DIAGRAMS_ROOT/2004_isabel_3D/diagrams.vtm" >/tmp/skot_matrix_verify.$$ 2>&1
  else
    env -u TTK_PLUGIN "$pv" "$MODULE_DIR/verify_ttk_matrices.py" \
      --input "$DIAGRAMS_ROOT/2004_isabel_3D/diagrams.vtm" >/tmp/skot_matrix_verify.$$ 2>&1
  fi
}

if ! verify_candidate "$PVPYTHON_PATH" "$TTK_PLUGIN_PATH"; then
  found=0
  mapfile -t PV_CANDIDATES < <(
    {
      command -v pvpython 2>/dev/null || true
      printf '%s\n' /usr/local/bin/pvpython /usr/bin/pvpython
      find "$HOME" /opt /usr/local -maxdepth 9 -type f -name pvpython -perm -u+x 2>/dev/null || true
    } | awk 'NF && !seen[$0]++'
  )
  mapfile -t PLUGIN_CANDIDATES < <(
    {
      printf '%s\n' \
        /usr/local/bin/plugins/TopologyToolKit/TopologyToolKit.so \
        /usr/local/lib/paraview*/plugins/TopologyToolKit/TopologyToolKit.so
      find "$HOME" /opt /usr/local -maxdepth 11 -type f \
        \( -name TopologyToolKit.so -o -name libTopologyToolKit.so \) 2>/dev/null || true
    } | awk 'NF && !seen[$0]++'
  )
  for pv in "${PV_CANDIDATES[@]}"; do
    if verify_candidate "$pv" ""; then
      PVPYTHON_PATH="$pv"; TTK_PLUGIN_PATH=""; found=1; break
    fi
    for plugin in "${PLUGIN_CANDIDATES[@]}"; do
      if verify_candidate "$pv" "$plugin"; then
        PVPYTHON_PATH="$pv"; TTK_PLUGIN_PATH="$plugin"; found=1; break 2
      fi
    done
  done
  rm -f /tmp/skot_matrix_verify.$$
  [ "$found" -eq 1 ] || fail "Could not find a pvpython executable loading the custom TTK filter."
fi
rm -f /tmp/skot_matrix_verify.$$

{
  printf 'PVPYTHON_PATH=%q\n' "$PVPYTHON_PATH"
  printf 'TTK_PLUGIN_PATH=%q\n' "$TTK_PLUGIN_PATH"
} > "$CONFIG_FILE"

export TTK_PLUGIN="$TTK_PLUGIN_PATH"
echo
printf 'MATRIX COMPUTATION — 12 COLLECTIONS\n'
printf 'pvpython : %s\n' "$PVPYTHON_PATH"
printf 'plugin   : %s\n' "${TTK_PLUGIN_PATH:-loaded automatically}"
if [ "$TIMING_RUNS" -eq 1 ]; then
  printf 'timings  : one run for d_SK,L=30, W_Gamma,L=30, and W2\n'
else
  printf 'timings  : mean over %s runs for d_SK,L=30, W_Gamma,L=30, and W2\n' "$TIMING_RUNS"
fi
printf 'DeltaLim : %s\n' "$DELTA_LIM"
printf 'output   : %s\n\n' "$OUTPUT_ROOT"

python3 -u "$MODULE_DIR/orchestrate_matrices.py" \
  --pvpython "$PVPYTHON_PATH" \
  --plugin "$TTK_PLUGIN_PATH" \
  --diagrams-root "$DIAGRAMS_ROOT" \
  --output-root "$OUTPUT_ROOT" \
  --compute-script "$MODULE_DIR/compute_collection_matrices.py" \
  --timing-runs "$TIMING_RUNS" \
  --delta-lim "$DELTA_LIM" \
  --heartbeat-seconds "$HEARTBEAT_SECONDS" \
  --collection-timeout-minutes "$COLLECTION_TIMEOUT_MINUTES"

python3 -u "$MODULE_DIR/finalize_matrices.py" \
  --output-root "$OUTPUT_ROOT" \
  --archive "$ARCHIVE" \
  --timing-runs "$TIMING_RUNS" \
  --delta-lim "$DELTA_LIM" \
  --pvpython "$PVPYTHON_PATH" \
  --plugin "$TTK_PLUGIN_PATH" \
  --diagrams-root "$DIAGRAMS_ROOT"

echo
echo "ALL TASKS COMPLETED."
echo "Matrix directory: $OUTPUT_ROOT"
echo "Archive for Stage 2: $ARCHIVE"
echo "Log: $LOG_FILE"
trap - ERR
