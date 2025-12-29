
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
  Rscript scripts/kinases/bin/generate_kinases_from_biomart.R --species human --out genesets/curated/kinases/outputs/kinases_human_union.csv
  ```
- Validate/export (GMT):
  ```bash
  Rscript scripts/kinases/bin/fetch_kinhub_and_merge.R
  Rscript scripts/kinases/bin/validations_index.R
  ```


## Checksums and temp outputs

- MD5 checksum files are stored alongside their parent CSVs (e.g. `genesets/curated/kinases/outputs/kinases_human_union.csv.md5`).
- Intermediate/temp CSVs are placed under `output/temp/` by default. Their accompanying `.md5` files are kept in the same folder (e.g. `output/temp/kinases_human_domain_interpro.csv.md5`).

Example filesystem layout:

```
genesets/curated/kinases/outputs/kinases_human_union.csv
genesets/curated/kinases/outputs/kinases_human_union.csv.md5
output/temp/kinases_human_domain_interpro.csv
output/temp/kinases_human_domain_interpro.csv.md5
output/temp/kinases_human_go_filtered.csv
output/temp/kinases_human_go_filtered.csv.md5
```

When regenerating outputs, write the checksum immediately next to the CSV to simplify discovery and manifest generation.

## Discovery / Indexing

- Use `scripts/kinases/bin/fetchers_index.R`, `scripts/kinases/bin/generate_kinases_from_biomart.R`, and `scripts/kinases/bin/validations_index.R` as short, discoverable entrypoints. They source helpers in `lib/`.

## Rationale and guidelines

- Keep `bin/` scripts thin (orchestration only). Move any shared logic into `lib/`.
- Archive originals rather than delete them immediately — they are available under `scripts/kinases/archive/`.
- After you confirm the canonical layout, we can remove archive contents or move them to an external snapshot.

Maintained by: bRNA3F AI Agent
