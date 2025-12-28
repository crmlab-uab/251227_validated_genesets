# Source Annotation Workflow (2025-12-27)

## Overview
This directory provides a fully reproducible workflow for generating authoritative, annotated source (kinase) tables for mouse and human. All intermediate scripts and files have been archived for clarity. Only the final, unified script and outputs remain.

## Usage
- Main script: `build_sources_annotation.R` (preferred) — a generalized builder that produces `kinases_human.csv` and `kinases_mouse.csv`.
- Run for mouse:
  ```
  Rscript kinases/build_sources_annotation.R --species mouse
  ```
- Run for human:
  ```
  Rscript kinases/build_sources_annotation.R --species human
  ```
- Outputs:
  - `kinases_mouse.csv` — Mouse source list with group, metabolic, and lipid annotations
  - `kinases_human.csv` — Human source list with group, metabolic, and lipid annotations

## Features
- Fetches source gene lists from Ensembl/BioMart
- Annotates kinase group/family (HGNC, human only)
- Annotates metabolic and lipid kinases using KEGG and org.*.eg.db
- All code and outputs are version-controlled and archived

- Supports additional validation sources via `kinases/val_sources/` (CSV, GMT, HTML).
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
