# Author: C. Ryan Miller
# Created: 2025-12-28 02:55 CST
# Commit: d11b4416b9252abe6807758f54f1659f18675324

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
