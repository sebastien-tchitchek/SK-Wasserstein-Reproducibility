#!/usr/bin/env bash
set -Eeuo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$MODULE_DIR/../.." && pwd)"
SCRIPTS_DIR="$MODULE_DIR"
CACHE_DIR="$PROJECT_ROOT/.cache_wmt_reproduction"
REUSED_CACHE=""

ARCHIVE="$CACHE_DIR/wassersteinMergeTreesData.tar.xz"
DATA_DIR="$CACHE_DIR/extracted_data"
WORK_DIR="$CACHE_DIR/work"
OUTPUT_DIR="$PROJECT_ROOT/NORMALIZED_DIAGRAM_RESULTS"
LOG_FILE="$PROJECT_ROOT/01_diagrams.log"
CONFIG_FILE="$PROJECT_ROOT/ttk_configuration.sh"
FILE_ID="14dx5Fjw0PnBgYwzVkU7iD8jf5lW8ZNDf"
DRIVE_URL="https://drive.google.com/file/d/${FILE_ID}/view?usp=sharing"
VIDAL_FILE_ID="1fxPAUmYT4UAUHh7K0w83znckyw_NTLan"
VIDAL_DRIVE_URL="https://drive.google.com/file/d/${VIDAL_FILE_ID}/view?usp=sharing"
VIDAL_ARCHIVE="$DATA_DIR/vidal_supplementary_package"
MAIN_INVENTORY_FILE="$OUTPUT_DIR/MAIN_ARCHIVE_INVENTORY.txt"
VIDAL_INVENTORY_FILE="$OUTPUT_DIR/VIDAL_SUPPLEMENT_INVENTORY.txt"
LOW_DISK_MAIN_REPORT="$OUTPUT_DIR/MAIN_LOW_DISK_REPORT.txt"
LOW_DISK_VIDAL_REPORT="$OUTPUT_DIR/VIDAL_LOW_DISK_REPORT.txt"

mkdir -p "$CACHE_DIR" "$WORK_DIR"
touch "$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

if [ -t 1 ]; then
  BOLD='\033[1m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'; RED='\033[0;31m'; RESET='\033[0m'
else
  BOLD=''; GREEN=''; YELLOW=''; RED=''; RESET=''
fi

say() { printf '%b\n' "$*"; }
fail() {
  say "${RED}${BOLD}ERROR:${RESET} $*"
  say "Details are recorded in: $LOG_FILE"
  exit 1
}

on_error() {
  local code=$?
  say ""
  say "${RED}${BOLD}The script stopped before completion.${RESET}"
  say "Rerun exactly: bash run.sh"
  say "Completed collections will be preserved."
  say "Log: $LOG_FILE"
  exit "$code"
}
trap on_error ERR

say ""
say "${BOLD}COMPUTATION OF NORMALIZED PERSISTENCE DIAGRAMS FOR SK WASSERSTEIN${RESET}"
say "- automatic download of the official scalar fields"
say "- classical TTK diagrams, all dimensions, DMS backend"
say "- one common normalization per dataset in [0,1]^2"
say "- no distance matrix is computed at this stage"
say "- processing one dataset at a time"
say "- automatic deletion of scalar fields after each dataset is validated"
say "- the large temporal dataset is opened only after disk space has been released"
say "- download of the Vidal supplement only when missing"
if [ -n "$REUSED_CACHE" ]; then
  say "${GREEN}Previously downloaded data were detected automatically.${RESET}"
  say "  Reused cache: $REUSED_CACHE"
fi
say ""

all_results_complete() {
  [ -f "$OUTPUT_DIR/ALL_DONE.txt" ] || return 1
  local markers manifests diagrams
  markers="$(find "$OUTPUT_DIR" -mindepth 2 -maxdepth 2 -type f -name 'DONE.ok' 2>/dev/null | wc -l)"
  manifests="$(find "$OUTPUT_DIR" -mindepth 2 -maxdepth 2 -type f -name 'diagrams.vtm' 2>/dev/null | wc -l)"
  diagrams="$(find "$OUTPUT_DIR" -mindepth 2 -maxdepth 2 -type f -name '*.vtu' 2>/dev/null | wc -l)"
  [ "$markers" -eq 12 ] && [ "$manifests" -eq 12 ] && [ "$diagrams" -eq 227 ]
}

if all_results_complete; then
  say "${GREEN}${BOLD}All 227 persistence diagrams are already present.${RESET}"
  say "Results: $OUTPUT_DIR"
  say "Nothing will be downloaded or recomputed."
  trap - ERR
  exit 0
fi

install_if_missing() {
  local missing=()
  command -v python3 >/dev/null 2>&1 || missing+=(python3)
  command -v curl >/dev/null 2>&1 || missing+=(curl)
  command -v tar >/dev/null 2>&1 || missing+=(tar)
  command -v xz >/dev/null 2>&1 || missing+=(xz-utils)
  if [ ${#missing[@]} -eq 0 ]; then
    return
  fi
  say "Automatically installing missing small utilities: ${missing[*]}"
  command -v sudo >/dev/null 2>&1 || fail "sudo is required to install ${missing[*]}."
  sudo apt-get update
  sudo apt-get install -y python3 curl tar xz-utils ca-certificates
}
install_if_missing

main_results_complete() {
  local key expected count
  while IFS=':' read -r key expected; do
    [ -f "$OUTPUT_DIR/$key/DONE.ok" ] || return 1
    [ -f "$OUTPUT_DIR/$key/diagrams.vtm" ] || return 1
    count="$(find "$OUTPUT_DIR/$key" -maxdepth 1 -type f -name '*.vtu' 2>/dev/null | wc -l)"
    [ "$count" -eq "$expected" ] || return 1
  done <<'EOF_MAIN_DATASETS'
2006_earthquake_3D:12
2008_ionization_front_2D:16
2008_ionization_front_3D:16
2014_volcanic_eruptions_2D:12
2016_viscous_fingering_3D:15
2017_cloud_processes_2D:12
2018_asteroid_impact_3D_clustering:7
2018_asteroid_impact_3D_temporal_subsampling:20
EOF_MAIN_DATASETS
  return 0
}

data_is_ready() {
  [ -f "$DATA_DIR/.extraction_complete" ] || return 1
  [ -f "$DATA_DIR/.main_collections_complete" ] && return 0
  [ -n "$(find "$DATA_DIR" -mindepth 1 -type f \
                 ! -name '.extraction_complete' -print -quit 2>/dev/null)" ]
}

archive_is_valid() {
  if data_is_ready && [ -s "$ARCHIVE" ]; then
    return 0
  fi
  [ -s "$ARCHIVE" ] && tar -tJf "$ARCHIVE" >/dev/null 2>&1
}

fallback_gdown() {
  local file_id="$1"
  local destination="$2"
  local drive_url="$3"
  local venv="$CACHE_DIR/venv_gdown"
  say "Fallback download method (gdown)..."
  if ! python3 -m venv "$venv" >/dev/null 2>&1; then
    say "Installing python3-venv..."
    sudo apt-get update
    sudo apt-get install -y python3-venv
    python3 -m venv "$venv"
  fi
  "$venv/bin/python" -m pip install --upgrade pip gdown
  rm -f "$destination" "$destination.part"
  "$venv/bin/gdown" --fuzzy "$drive_url" -O "$destination"
}

download_archive() {
  if data_is_ready || main_results_complete; then
    say "${GREEN}Existing data will be reused: no new download.${RESET}"
    return
  fi
  if archive_is_valid; then
    say "${GREEN}The official archive is already present and valid.${RESET}"
    return
  fi

  rm -f "$ARCHIVE" "$ARCHIVE.part"
  say "Downloading the official WassersteinMergeTreesData archive..."
  local direct_url="https://drive.usercontent.google.com/download?id=${FILE_ID}&export=download&confirm=t"
  if curl -L --fail --retry 6 --retry-delay 3 --retry-all-errors \
      --connect-timeout 30 --output "$ARCHIVE.part" "$direct_url"; then
    mv "$ARCHIVE.part" "$ARCHIVE"
  fi

  if ! archive_is_valid; then
    rm -f "$ARCHIVE" "$ARCHIVE.part"
    say "The direct method was insufficient; using the included Google Drive downloader..."
    python3 "$SCRIPTS_DIR/download_google_drive.py" "$FILE_ID" "$ARCHIVE" || true
  fi

  if ! archive_is_valid; then
    fallback_gdown "$FILE_ID" "$ARCHIVE" "$DRIVE_URL"
  fi

  archive_is_valid || fail "The downloaded archive is not a valid tar.xz file."
  say "${GREEN}Download completed and verified.${RESET}"
}

extract_archive() {
  if xz --robot --list "$ARCHIVE" >/tmp/wmt_xz_info.$$ 2>/dev/null; then
    local uncompressed_bytes free_kb needed_kb
    uncompressed_bytes="$(awk -F '\t' '$1=="totals" {print $5}' /tmp/wmt_xz_info.$$ | tail -1)"
    rm -f /tmp/wmt_xz_info.$$
    if [[ "$uncompressed_bytes" =~ ^[0-9]+$ ]]; then
      free_kb="$(df -Pk "$CACHE_DIR" | awk 'NR==2 {print $4}')"
      needed_kb=$(( uncompressed_bytes / 1024 + 2 * 1024 * 1024 ))
      if [ "$free_kb" -lt "$needed_kb" ]; then
        fail "Insufficient disk space to extract the data. Free some space and rerun the script."
      fi
    fi
  fi

  if ! data_is_ready && ! main_results_complete; then
    say "Extracting the first data layer..."
    rm -rf "$DATA_DIR"
    mkdir -p "$DATA_DIR"
    tar -xJf "$ARCHIVE" -C "$DATA_DIR"
    touch "$DATA_DIR/.extraction_complete"
    say "${GREEN}First layer extracted.${RESET}"
  else
    mkdir -p "$DATA_DIR"
    [ -f "$DATA_DIR/.extraction_complete" ] || touch "$DATA_DIR/.extraction_complete"
    say "${GREEN}Existing data found and reused.${RESET}"
  fi
}



recognized_archive_is_valid() {
  local path="$1"
  [ -s "$path" ] || return 1
  python3 - "$SCRIPTS_DIR" "$path" <<'PYARCHIVE'
from pathlib import Path
import sys
sys.path.insert(0, sys.argv[1])
from prepare_archives import archive_kind
sys.exit(0 if archive_kind(Path(sys.argv[2])) is not None else 1)
PYARCHIVE
}

vidal_files_present() {
  local name
  for name in \
    isabella_velocity_goodEnsemble.vti \
    seaSurfaceHeightGoodEnsemble.vti \
    startingVortexGoodEnsemble.vti \
    vortexStreetGoodEnsemble2.vti; do
    [ -n "$(find "$DATA_DIR" -type f -iname "$name" -print -quit 2>/dev/null)" ] || return 1
  done
  return 0
}

prepare_nested_archives() {
  mkdir -p "$OUTPUT_DIR"
  say ""
  say "Opening the second archive layer from the previously downloaded large file..."
  python3 "$SCRIPTS_DIR/prepare_archives.py" \
    --data-root "$DATA_DIR" \
    --outer-archive "$ARCHIVE" \
    --report "$MAIN_INVENTORY_FILE"
}

download_vidal_archive() {
  if recognized_archive_is_valid "$VIDAL_ARCHIVE"; then
    say "The Isabel/Vortex/SSH supplementary package is already present and valid."
    return
  fi

  rm -f "$VIDAL_ARCHIVE" "$VIDAL_ARCHIVE.part"
  say "Downloading the official Isabel / Vortex / SSH supplement..."
  local direct_url="https://drive.usercontent.google.com/download?id=${VIDAL_FILE_ID}&export=download&confirm=t"
  if curl -L --fail --retry 6 --retry-delay 3 --retry-all-errors \
      --connect-timeout 30 --output "$VIDAL_ARCHIVE.part" "$direct_url"; then
    mv "$VIDAL_ARCHIVE.part" "$VIDAL_ARCHIVE"
  fi

  if ! recognized_archive_is_valid "$VIDAL_ARCHIVE"; then
    rm -f "$VIDAL_ARCHIVE" "$VIDAL_ARCHIVE.part"
    say "The direct method was insufficient; using the included Google Drive downloader..."
    python3 "$SCRIPTS_DIR/download_google_drive.py" \
      "$VIDAL_FILE_ID" "$VIDAL_ARCHIVE" || true
  fi

  if ! recognized_archive_is_valid "$VIDAL_ARCHIVE"; then
    fallback_gdown "$VIDAL_FILE_ID" "$VIDAL_ARCHIVE" "$VIDAL_DRIVE_URL"
  fi

  recognized_archive_is_valid "$VIDAL_ARCHIVE" \
    || fail "The downloaded Isabel/Vortex/SSH supplement is not a recognized archive."
  say "Supplementary download completed and verified."
}

ensure_vidal_data() {
  if vidal_files_present; then
    say "The four Isabel/Vortex/SSH collections are already present."
    return
  fi

  download_vidal_archive
  say "Opening the supplementary package..."
  python3 "$SCRIPTS_DIR/prepare_archives.py" \
    --data-root "$DATA_DIR" \
    --report "$VIDAL_INVENTORY_FILE"

  vidal_files_present \
    || fail "The four Isabel/Vortex/SSH files are still missing. See $VIDAL_INVENTORY_FILE."
  say "The four supplementary collections have been prepared."
}

verify_prepared_data() {
  local count
  count="$(find "$DATA_DIR" -type f \
    \( -iname '*.vtu' -o -iname '*.vti' -o -iname '*.vtk' -o -iname '*.pvtu' -o -iname '*.pvti' \) \
    2>/dev/null | wc -l)"
  if [ "$count" -eq 0 ]; then
    fail "No VTK field was found after extraction. See $MAIN_INVENTORY_FILE."
  fi
  say "Prepared data: $count VTK files detected before selection of the 227 members."
}

if [ -f "$CONFIG_FILE" ]; then
  source "$CONFIG_FILE"
fi
PVPYTHON_PATH="${PVPYTHON_PATH:-${PVPYTHON:-}}"
TTK_PLUGIN_PATH="${TTK_PLUGIN_PATH:-${TTK_PLUGIN:-}}"

unique_append() {
  local value="$1"; shift
  [ -n "$value" ] || return 0
  [ -x "$value" ] || return 0
  local existing
  for existing in "${PVPYTHON_CANDIDATES[@]:-}"; do
    [ "$existing" = "$value" ] && return 0
  done
  PVPYTHON_CANDIDATES+=("$value")
}

unique_plugin_append() {
  local value="$1"; shift
  [ -n "$value" ] || return 0
  [ -f "$value" ] || return 0
  local existing
  for existing in "${PLUGIN_CANDIDATES[@]:-}"; do
    [ "$existing" = "$value" ] && return 0
  done
  PLUGIN_CANDIDATES+=("$value")
}

find_plugins() {
  PLUGIN_CANDIDATES=()
  unique_plugin_append "$TTK_PLUGIN_PATH"

  local plugin_dir plugin
  IFS=':' read -r -a _pv_plugin_dirs <<< "${PV_PLUGIN_PATH:-}"
  for plugin_dir in "${_pv_plugin_dirs[@]:-}"; do
    [ -d "$plugin_dir" ] || continue
    while IFS= read -r plugin; do unique_plugin_append "$plugin"; done < <(
      find "$plugin_dir" -maxdepth 5 -type f \
        \( -name 'TopologyToolKit.so' -o -name 'libTopologyToolKit.so' \) \
        2>/dev/null | sort -u
    )
  done

  while IFS= read -r plugin; do unique_plugin_append "$plugin"; done < <(
    find "$HOME" /opt /usr/local -maxdepth 11 -type f \
      \( -name 'TopologyToolKit.so' -o -name 'libTopologyToolKit.so' \) \
      -not -path '*/.cache/*' -not -path '*/Trash/*' 2>/dev/null | sort -u
  )
}

check_pvpython() {
  local candidate="$1"
  local plugin="${2:-}"
  if [ -n "$plugin" ]; then
    TTK_PLUGIN="$plugin" "$candidate" "$SCRIPTS_DIR/verify_ttk.py" >/tmp/wmt_ttk_check.$$ 2>&1
  else
    env -u TTK_PLUGIN "$candidate" "$SCRIPTS_DIR/verify_ttk.py" >/tmp/wmt_ttk_check.$$ 2>&1
  fi
}

find_pvpython() {
  PVPYTHON_CANDIDATES=()
  unique_append "$PVPYTHON_PATH"
  local command_name located
  for command_name in pvpython pvpython5.13 pvpython5.12 pvpython5.11 pvpython5.10; do
    located="$(command -v "$command_name" 2>/dev/null || true)"
    unique_append "$located"
  done

  while IFS= read -r located; do unique_append "$located"; done < <(
    find "$HOME" /opt /usr/local -maxdepth 9 -type f -name pvpython -perm -u+x \
      -not -path '*/.cache/*' -not -path '*/Trash/*' 2>/dev/null | sort -u
  )
  find_plugins

  local candidate plugin
  local fallback_pv=""
  local fallback_plugin=""
  for candidate in "${PVPYTHON_CANDIDATES[@]}"; do
    if check_pvpython "$candidate" ""; then
      if grep -q '^TTK_SKOT_CUSTOM=1' /tmp/wmt_ttk_check.$$; then
        PVPYTHON_PATH="$candidate"
        TTK_PLUGIN_PATH=""
        rm -f /tmp/wmt_ttk_check.$$
        return 0
      fi
      if [ -z "$fallback_pv" ]; then
        fallback_pv="$candidate"
        fallback_plugin=""
      fi
    fi

    for plugin in "${PLUGIN_CANDIDATES[@]:-}"; do
      if check_pvpython "$candidate" "$plugin"; then
        if grep -q '^TTK_SKOT_CUSTOM=1' /tmp/wmt_ttk_check.$$; then
          PVPYTHON_PATH="$candidate"
          TTK_PLUGIN_PATH="$plugin"
          rm -f /tmp/wmt_ttk_check.$$
          return 0
        fi
        if [ -z "$fallback_pv" ]; then
          fallback_pv="$candidate"
          fallback_plugin="$plugin"
        fi
      fi
    done
  done
  rm -f /tmp/wmt_ttk_check.$$
  if [ -n "$fallback_pv" ]; then
    PVPYTHON_PATH="$fallback_pv"
    TTK_PLUGIN_PATH="$fallback_plugin"
    return 0
  fi
  return 1
}

if ! find_pvpython; then
  say ""
  say "${YELLOW}The ParaView/TTK pvpython executable could not be detected automatically.${RESET}"
  say "The TTK source ZIP is not executable by itself: the custom version must already be compiled."
  if [ ! -t 0 ]; then
    fail "Rerun in a terminal and provide the requested path to pvpython."
  fi
  read -r -p "Enter the full path to the pvpython executable, then press Enter: " PVPYTHON_PATH
  PVPYTHON_PATH="${PVPYTHON_PATH/#\~/$HOME}"
  [ -x "$PVPYTHON_PATH" ] || fail "This pvpython file does not exist or is not executable."

  if ! check_pvpython "$PVPYTHON_PATH" "$TTK_PLUGIN_PATH"; then
    say "The TTK filter is not loaded by this pvpython executable."
    read -r -p "Path to TopologyToolKit.so (or press Enter if none is needed): " TTK_PLUGIN_PATH
    TTK_PLUGIN_PATH="${TTK_PLUGIN_PATH/#\~/$HOME}"
    check_pvpython "$PVPYTHON_PATH" "$TTK_PLUGIN_PATH" || fail "This pvpython executable cannot load TTK PersistenceDiagram."
  fi
fi

say "selected pvpython: $PVPYTHON_PATH"
if [ -n "$TTK_PLUGIN_PATH" ]; then
  say "selected TTK plugin: $TTK_PLUGIN_PATH"
fi

if [ -n "$TTK_PLUGIN_PATH" ]; then
  TTK_PLUGIN="$TTK_PLUGIN_PATH" "$PVPYTHON_PATH" "$SCRIPTS_DIR/verify_ttk.py" --full
else
  env -u TTK_PLUGIN "$PVPYTHON_PATH" "$SCRIPTS_DIR/verify_ttk.py" --full
fi
say "${GREEN}TTK test succeeded.${RESET}"

{
  printf 'PVPYTHON_PATH=%q\n' "$PVPYTHON_PATH"
  printf 'TTK_PLUGIN_PATH=%q\n' "$TTK_PLUGIN_PATH"
} > "$CONFIG_FILE"

download_archive
extract_archive
find "$DATA_DIR" -type f \
  \( -iname '*2015_dark_matter_3D*.tar.xz' -o -iname '*2015_dark_matter_3D*.txz' \) \
  -print -delete 2>/dev/null || true
mkdir -p "$OUTPUT_DIR"

say ""
say "${BOLD}LOW-DISK MODE: main collections processed one at a time.${RESET}"
say "After each collection is validated, its scalar fields are deleted."
say "The normalized results are preserved."
say ""

python3 -u "$SCRIPTS_DIR/orchestrate_low_space.py" \
  --phase main \
  --data-root "$DATA_DIR" \
  --output-root "$OUTPUT_DIR" \
  --work-root "$WORK_DIR" \
  --pvpython "$PVPYTHON_PATH" \
  --plugin "$TTK_PLUGIN_PATH" \
  --compute-script "$SCRIPTS_DIR/compute_one_collection.py" \
  --outer-archive "$ARCHIVE" \
  --report "$LOW_DISK_MAIN_REPORT"

say ""
say "The eight main datasets are complete; preparing the Vidal supplement..."
ensure_vidal_data

python3 -u "$SCRIPTS_DIR/orchestrate_low_space.py" \
  --phase vidal \
  --data-root "$DATA_DIR" \
  --output-root "$OUTPUT_DIR" \
  --work-root "$WORK_DIR" \
  --pvpython "$PVPYTHON_PATH" \
  --plugin "$TTK_PLUGIN_PATH" \
  --compute-script "$SCRIPTS_DIR/compute_one_collection.py" \
  --report "$LOW_DISK_VIDAL_REPORT"

python3 -u "$SCRIPTS_DIR/verify_final_results.py" --output-root "$OUTPUT_DIR"

say ""
say "${GREEN}${BOLD}ALL TASKS COMPLETED.${RESET}"
say "The 227 normalized persistence diagrams are in:"
say "  $OUTPUT_DIR"
say ""
say "For each dataset, open its diagrams.vtm file in ParaView."
say "The matrices are computed automatically by the next stage of the main launcher."

say "Removing large temporary files that are no longer needed..."
rm -rf "$DATA_DIR" "$WORK_DIR" "$CACHE_DIR/venv_gdown"
rm -f "$ARCHIVE" "$ARCHIVE.part" "$VIDAL_ARCHIVE" "$VIDAL_ARCHIVE.part"
say "${GREEN}Cleanup completed.${RESET}"
say "Complete log: $LOG_FILE"
trap - ERR
