#!/usr/bin/env bash
# run_pipeline.sh
# Author: C. Ryan Miller
# Purpose: Run the kinases pipeline end-to-end
# Usage: ./run_pipeline.sh [--species=human|mouse] [--skip-fetch] [--skip-validate]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Defaults
SPECIES="human"
SKIP_FETCH=false
SKIP_VALIDATE=false
VALIDATE_ARGS=""

# Parse arguments
for arg in "$@"; do
  case $arg in
    --species=*)
      SPECIES="${arg#*=}"
      ;;
    --skip-fetch)
      SKIP_FETCH=true
      ;;
    --skip-validate)
      SKIP_VALIDATE=true
      ;;
    --skip-genome)
      VALIDATE_ARGS="$VALIDATE_ARGS --skip-genome"
      ;;
    --skip-comprehensive)
      VALIDATE_ARGS="$VALIDATE_ARGS --skip-comprehensive"
      ;;
    --help|-h)
      echo "Usage: $0 [--species=human|mouse] [--skip-fetch] [--skip-validate]"
      echo ""
      echo "Options:"
      echo "  --species=human|mouse   Species to process (default: human)"
      echo "  --skip-fetch            Skip BioMart fetch (use existing inputs)"
      echo "  --skip-validate         Skip validation step entirely"
      echo "  --skip-genome           Skip genome GTF validation (faster)"
      echo "  --skip-comprehensive    Skip comprehensive UniProt validation (faster)"
      echo ""
      exit 0
      ;;
  esac
done

cd "$REPO_ROOT"
echo "=== Kinases Pipeline ==="
echo "Species: $SPECIES"
echo "Repo root: $REPO_ROOT"
echo ""

# Step 01: Fetch from BioMart
if [ "$SKIP_FETCH" = false ]; then
  echo "[Step 01] Fetching kinases from BioMart..."
  Rscript scripts/kinases/01_fetch_geneset_BioMart.R --species="$SPECIES"
else
  echo "[Step 01] Skipped (--skip-fetch)"
fi

# Step 02: Reconcile Manning vs KinHub sources
echo ""
echo "[Step 02] Reconciling Manning vs KinHub sources..."
Rscript scripts/kinases/02_reconcile_sources.R

# Step 03: Build annotations
echo ""
echo "[Step 03] Building annotations..."
Rscript scripts/kinases/03_build_annotations.R --species="$SPECIES"

# Step 04: Map human to mouse (only for human)
if [ "$SPECIES" = "human" ]; then
  echo ""
  echo "[Step 04] Mapping human kinases to mouse orthologs..."
  Rscript scripts/kinases/04_map_human_to_mouse.R
fi

# Step 05: Export GMT files
echo ""
echo "[Step 05] Exporting GMT files..."
Rscript scripts/kinases/05_export_gmt.R

# Step 06: Validation (optional)
if [ "$SKIP_VALIDATE" = false ]; then
  echo ""
  echo "[Step 06] Running validation..."
  # shellcheck disable=SC2086
  Rscript scripts/kinases/06_validate.R $VALIDATE_ARGS
else
  echo ""
  echo "[Step 06] Skipped (--skip-validate)"
fi

echo ""
echo "=== Pipeline Complete ==="
echo "Outputs in: curated/kinases/outputs/"
ls -la curated/kinases/outputs/*.gmt curated/kinases/outputs/*.csv 2>/dev/null | head -10

echo ""
echo "=== CORAL Kinome Tree (Optional) ==="
echo "To generate CORAL visualization input from DESeq2 results:"
echo "  Rscript scripts/kinases/07_generate_coral_tree.R --log2fc=<DEG_all.csv> --comparison=<name>"
echo ""
echo "Example with 2x2 factorial Driver comparison:"
echo "  Rscript scripts/kinases/07_generate_coral_tree.R \\"
echo "    --log2fc=/data/251207_mNSC/output_factorial_2x2/csv/Driver_Pdgfra_vs_Egfr_Host_BL6_TSG_Y_Model_NS2_protein_coding_additive/251224_Driver_Pdgfra_vs_Egfr_Host_BL6_TSG_Y_Model_NS2_protein_coding_additive_DriverPdgfra_DEG_all.csv \\"
echo "    --comparison='Driver_Pdgfra_vs_Egfr' \\"
echo "    --output-prefix='mNSC_NS2'"
echo ""
echo "Then load output in CORAL: https://phanstiel-lab.shinyapps.io/CORAL/"
