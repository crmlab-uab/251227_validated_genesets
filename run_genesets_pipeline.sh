#!/usr/bin/env bash
set -euo pipefail

# Simple wrapper to run the R-based genesets pipeline
# Usage: ./run_genesets_pipeline.sh [--species human|mouse] [extra args]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found on PATH" >&2
  exit 2
fi

echo "Running genesets pipeline in $PWD"
Rscript run_genesets_pipeline.R "$@"
