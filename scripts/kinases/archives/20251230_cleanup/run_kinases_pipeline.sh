#!/usr/bin/env bash
set -euo pipefail

usage(){
  cat <<EOF
Usage: $(basename "$0") [--from N] [--to N] [--dry-run] [--force] [--log-dir DIR]

Runs kinases scripts in order: executes R scripts named `NN_*.R` in this folder.
Options:
  --from N     Start step index (1-based)
  --to N       End step index (inclusive)
  --dry-run    Print commands but don't execute
  --force      Continue even if outputs exist
  --log-dir DIR Write logs to DIR (default: sessions/<timestamp>_kinases_run)
EOF
  exit 1
}

FROM=1
TO=999
DRY_RUN=0
FORCE=0
LOG_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="$2"; shift 2;;
    --to) TO="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    --force) FORCE=1; shift;;
    --log-dir) LOG_DIR="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

if [[ -z "$LOG_DIR" ]]; then
  TS="$(date +"%y%m%d_%H%M%S")"
  LOG_DIR="$REPO_ROOT/sessions/${TS}_kinases_run"
fi
mkdir -p "$LOG_DIR"


# discover numbered R scripts in this bin dir and sort by numeric prefix
mapfile -t SCRIPTS < <(
  find "$SCRIPT_DIR" -maxdepth 1 -type f -name '[0-9][0-9]_*.R' -printf '%f\n' \
  | awk -F'_' '{print $1 "\t" $0}' \
  | sort -n -k1,1 \
  | cut -f2 \
  | sed "s|^|$SCRIPT_DIR/|"
)
if [[ ${#SCRIPTS[@]} -eq 0 ]]; then
  echo "No numbered R scripts found in $SCRIPT_DIR" | tee -a "$LOG_DIR/run.log"
  exit 1
fi

TOTAL=${#SCRIPTS[@]}
if (( FROM < 1 )); then FROM=1; fi
if (( TO > TOTAL )); then TO=$TOTAL; fi

echo "Repository: $REPO_ROOT" | tee -a "$LOG_DIR/run.log"
echo "Scripts dir: $SCRIPT_DIR" | tee -a "$LOG_DIR/run.log"
echo "Running steps $FROM..$TO of $TOTAL" | tee -a "$LOG_DIR/run.log"

STEP=0
for SCRIPT in "${SCRIPTS[@]}"; do
  (( STEP++ ))
  if (( STEP < FROM || STEP > TO )); then
    continue
  fi
  NAME="$(basename "$SCRIPT")"
  echo "---" | tee -a "$LOG_DIR/run.log"
  echo "Step $STEP/$TOTAL: $NAME" | tee -a "$LOG_DIR/run.log"
  echo "Working dir: $REPO_ROOT" | tee -a "$LOG_DIR/run.log"
  CMD=(Rscript "$SCRIPT")
  echo "Command: ${CMD[*]}" | tee -a "$LOG_DIR/run.log"

  if (( DRY_RUN )); then
    echo "Dry-run: skipping execution" | tee -a "$LOG_DIR/run.log"
    continue
  fi

  # run script from repository root to keep relative paths stable
  if ! (cd "$REPO_ROOT" && "${CMD[@]}" 2>&1 | tee -a "$LOG_DIR/${STEP}_${NAME}.log"); then
    echo "Step $STEP failed. See $LOG_DIR/${STEP}_${NAME}.log" | tee -a "$LOG_DIR/run.log"
    exit 1
  fi

  echo "Step $STEP completed" | tee -a "$LOG_DIR/run.log"
done

echo "All requested steps completed successfully." | tee -a "$LOG_DIR/run.log"
