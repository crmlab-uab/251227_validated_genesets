
# Kinases scripts

Canonical layout (current):

- `bin/` – executable entrypoints and thin runners (index/fetch/generate/validate)
- `lib/` – reusable library functions (mapping, cleaning, merging, helpers)
- `utils/` – small utilities and helpers
- `archive/` – deprecated or original scripts (moved here on 2025-12-28)
- `temp/` – runtime intermediate outputs (do not commit; keep checksums instead)

Use `bin/` scripts as run entrypoints; they should source `lib/` functions rather than duplicating logic.

## Status (2025-12-28)

- Original, top-level kinase scripts were moved into `scripts/kinases/archive/` to preserve history while making `bin/` and `lib/` the canonical layout. See `scripts/kinases/archive/` for originals.

## Recommended quick-start

- Fetch sources (example):
  ```bash
  Rscript scripts/kinases/bin/fetch_kinome.R --species human
  ```
- Generate canonical kinases list from BioMart (GO + InterPro union):
  ```bash
  Rscript scripts/kinases/bin/generate_kinases_from_biomart.R --species human --out genesets/curated/kinases/kinases_human_union.csv
  ```
- Validate/export (GMT):
  ```bash
  Rscript scripts/kinases/bin/fetch_kinhub_and_merge.R
  Rscript scripts/kinases/bin/validations_index.R
  ```

## Checksums and manifest

- Recompute checksums for canonical outputs (example):
  ```bash
  cd genesets/curated/kinases
  mkdir -p checksums
  md5sum kinases_human_union.csv > checksums/kinases_human_union.csv.md5
  md5sum temp/kinases_human_domain.csv > temp/checksums/kinases_human_domain.csv.md5
  md5sum temp/kinases_human_go.csv > temp/checksums/kinases_human_go.csv.md5
  ```

## Discovery / Indexing

- Use `scripts/kinases/bin/fetchers_index.R`, `scripts/kinases/bin/generate_kinases_from_biomart.R`, and `scripts/kinases/bin/validations_index.R` as short, discoverable entrypoints. They source helpers in `lib/`.

## Rationale and guidelines

- Keep `bin/` scripts thin (orchestration only). Move any shared logic into `lib/`.
- Archive originals rather than delete them immediately — they are available under `scripts/kinases/archive/`.
- After you confirm the canonical layout, we can remove archive contents or move them to an external snapshot.

Maintained by: bRNA3F AI Agent
