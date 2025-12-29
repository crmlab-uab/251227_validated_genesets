# Source Annotation Workflow (2025-12-27)

## Overview
This directory provides a fully reproducible workflow for generating authoritative, annotated source (kinase) tables for mouse and human. All intermediate scripts and files have been archived for clarity. Only the final, unified script and outputs remain.

## Usage
- Recommended entry points (grouped index scripts):
  - `Rscript scripts/kinases/fetchers_index.R` — source all fetcher scripts (downloads and fetch wrappers)
  - `Rscript scripts/kinases/mappings_index.R` — source mapping/merge helpers (uniprot/biomart mappings)
  - `Rscript scripts/kinases/annotations_index.R` — source annotation helpers (KEGG/Manning etc.)
  - `Rscript scripts/kinases/validations_index.R` — source validation and export helpers (validation, export GMT)

- Legacy main builder: `build_sources_annotation.R` remains and can be used to run the full pipeline; the index scripts provide a simpler way to load grouped helpers for interactive use.

- Example: run the full builder for mouse (legacy main script):
  ```
  Rscript scripts/kinases/build_sources_annotation.R --species mouse
  ```
- Outputs:
  - `kinases_mouse.csv` — Mouse source list with group, metabolic, and lipid annotations
  - `kinases_human.csv` — Human source list with group, metabolic, and lipid annotations

## Features
- Fetches source gene lists from Ensembl/BioMart
- Annotates kinase group/family (HGNC, human only)
- Annotates metabolic and lipid kinases using KEGG and org.*.eg.db
- All code and outputs are version-controlled and archived

- Supports additional validation sources via `genesets/curated/` (preferred). Legacy paths `val_sources/` and `kinases/val_sources/` are still recognized.
  - CSVs are merged by Ensembl ID or gene symbol when possible.
  - GMTs are parsed; gene sets are attached as a `val_sources` annotation.
  - HTML pages containing KinHub data will be delegated to the KinHub parser when the filename contains `kinhub`; otherwise the first HTML table is attempted.

## Directory Structure
- `build_sources_annotation.R` — Main script (preferred)
- `kinases_mouse.csv`, `kinases_human.csv` — Final outputs
 - Manning supplement CSV: `kinases/data/manning_2002_TableS1.csv` (moved from repo root)
- `archives/` — All legacy scripts and intermediate files
- `docs/` — Documentation (if present)

## Notes
- For GSEA or bRNA3F, use the output CSVs to generate GMT files as needed, ensuring gene symbols match your dataset species.
- For any updates, edit only the main script and re-run for both species.

---
Maintained by: bRNA3F AI Agent
