#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "$ROOT/.." && pwd)"
STAGE1="$PACKAGE_ROOT/01_COMPUTE_DIAGRAMS_AND_MATRICES"
OUTPUT="$PACKAGE_ROOT/results/analysis"
VENV="$ROOT/.venv_reproduction"
LOG="$ROOT/analysis.log"
EXTRACTED="$ROOT/.extracted_matrices"
touch "$LOG"
exec > >(tee -a "$LOG") 2>&1
fail() { echo "ERROR: $*" >&2; exit 1; }
MATRICES_DIR="${MATRICES_DIR:-}"
if [ -n "$MATRICES_DIR" ]; then
  MATRICES_DIR="$(realpath -m "$MATRICES_DIR")"
elif [ -d "$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS" ]; then
  MATRICES_DIR="$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS"
elif [ -f "$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS.zip" ]; then
  rm -rf "$EXTRACTED"
  mkdir -p "$EXTRACTED"
  python3 - "$STAGE1/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS.zip" "$EXTRACTED" <<'PY_EXTRACT'
import sys, zipfile
from pathlib import Path
archive, destination = map(Path, sys.argv[1:])
with zipfile.ZipFile(archive) as zf:
    zf.extractall(destination)
PY_EXTRACT
  MATRICES_DIR="$EXTRACTED/SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS"
else
  fail "No matrix was found."
fi
[ -d "$MATRICES_DIR" ] || fail "Matrix directory is missing: $MATRICES_DIR"
[ -f "$MATRICES_DIR/2004_isabel_3D/W2/matrix_distance.npy" ] || fail "Matrix directory is incomplete: $MATRICES_DIR"
command -v python3 >/dev/null 2>&1 || fail "python3 is required."
if python3 - <<'PYVERS' >/dev/null 2>&1
import numpy, pandas, scipy, sklearn, matplotlib
expected = ('2.3.5','2.2.3','1.17.0','1.8.0','3.10.8')
actual = (numpy.__version__, pandas.__version__, scipy.__version__, sklearn.__version__, matplotlib.__version__)
raise SystemExit(0 if actual == expected else 1)
PYVERS
then
  PYTHON_CMD="python3"
else
  if [ ! -d "$VENV" ]; then
    if ! python3 -m venv "$VENV"; then
      command -v sudo >/dev/null 2>&1 || fail "python3-venv is required."
      sudo apt-get update
      sudo apt-get install -y python3-venv
      python3 -m venv "$VENV"
    fi
  fi
  export PIP_DISABLE_PIP_VERSION_CHECK=1
  REQ_HASH="$(sha256sum "$ROOT/requirements.txt" | awk '{print $1}')"
  INSTALLED_HASH=""
  [ -f "$VENV/.requirements.sha256" ] && INSTALLED_HASH="$(cat "$VENV/.requirements.sha256")"
  if [ "$REQ_HASH" != "$INSTALLED_HASH" ]; then
    "$VENV/bin/python" -m pip install --only-binary=:all: -r "$ROOT/requirements.txt"
    printf '%s\n' "$REQ_HASH" > "$VENV/.requirements.sha256"
  fi
  PYTHON_CMD="$VENV/bin/python"
fi
rm -rf "$OUTPUT"
mkdir -p "$OUTPUT"
"$PYTHON_CMD" "$ROOT/scripts/reproduce_results.py" --matrices "$MATRICES_DIR" --output "$OUTPUT"
