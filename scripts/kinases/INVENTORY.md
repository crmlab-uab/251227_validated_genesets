# Kinases scripts inventory

This file lists each script under `scripts/kinases/` with a short purpose, typical inputs/outputs, and a recommended canonical destination within `scripts/kinases/`.

Format: Script — Purpose — Inputs — Outputs — Recommended destination

- `annotate_lipid_kinases_from_kegg.R` — Annotate kinases with lipid pathway membership via KEGG — inputs: `human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv` — outputs: same CSV updated — `lib`
- `annotate_metabolic_kinases_from_kegg.R` — Annotate kinases with metabolic pathway membership via KEGG — inputs: `human_kinome_with_group_family_GOlast_custom2.csv` — outputs: `human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv` — `lib`
- `annotations_index.R` — Index that sources annotation scripts — inputs: none (sources other scripts) — outputs: none — `bin`
- `augment_matching_with_aliases.R` — Augment baseline kinome with Manning aliases/HGNC lookups and write match CSV — inputs: `kinases_human.csv`, Manning table — outputs: `kinases_human.with_aliases_matches.csv` — `lib`
- `build_kinome_annotation.R` — Orchestrates BioMart fetch + HGNC group merge to produce human/mouse kinome CSVs — inputs: Manning file, BioMart — outputs: `kinases_human.csv`, `kinases_mouse.csv` — `bin`
- `clean_comprehensive_mouse_kinome.R` — Clean and enrich comprehensive mouse kinome export from BioMart — inputs: `comprehensive_mouse_kinome_biomart.csv` — outputs: `comprehensive_mouse_kinome_cleaned.csv`, `comprehensive_mouse_kinome_unique.csv` — `lib`
- `comprehensive_kinase_validation.R` — Validate mouse kinases across UniProt/GO/org.* DBs; produce validation CSV and report — inputs: curated kinases CSV — outputs: `mouse_kinome_validation_results.csv`, `VALIDATION_REPORT.md` — `bin`
- `export_kinase_gmt.R` — Export kinome sets as GMT files for GSEA — inputs: curated combined CSV — outputs: multiple `.gmt` files — `lib`
- `fetch_comprehensive_mouse_kinome_biomart.R` — BioMart query for mouse kinases by GO — inputs: BioMart — outputs: `comprehensive_mouse_kinome_biomart.csv` — `bin`
- `fetch_definitive_mouse_kinome.R` — Map curated human kinases to mouse via BioMart (definitive list) — inputs: `201006_composite_kinases_curated.csv` — outputs: `mouse_kinome_definitive.csv` — `bin`
- `fetch_hgnc_kinase_groups.R` — Fetch HGNC kinase group annotations via REST — inputs: network/HGNC REST — outputs: `hgnc_kinase_groups.csv` — `bin`
- `fetch_kinhub_and_merge.R` — Parse KinHub mapping and merge into baseline kinases CSV — inputs: `kinhub_mapping_raw.tsv`, baseline kinases — outputs: `kinases_human.with_kinhub.csv` — `lib`
- `fetch_kinome.R` — Generic BioMart fetcher for human/mouse kinome (GO-based) — inputs: BioMart — outputs: `human_kinome.csv` / `mouse_kinome.csv` — `bin`
- `fetch_mouse_kinome_from_kinhub.R` — Scrape KinHub and map to mouse orthologs via BioMart — inputs: KinHub HTML/table — outputs: `mouse_kinome_from_kinhub.csv` — `bin`
- `fetch_uniprot_kinase_mapping_api.R` — Query UniProt API for mouse kinase mapping — inputs: curated kinases CSV — outputs: `uniprot_mouse_kinase_idmapping.tab` — `lib`
- `fetch_uniprot_mouse_mapping_api.R` — General UniProt mouse mapping fetcher (paginated) — inputs: UniProt API — outputs: `uniprot_mouse_idmapping_selected.tab` — `lib`
- `fetch_val_sources_and_merge.R` — Generic merge of validation sources (CSV/XLS/GMT/HTML) into baseline kinases table — inputs: files under `val_sources/` — outputs: `kinases/kinases_human.with_val_sources.csv` — `lib`
- `fetchers_index.R` — Index to source fetcher entrypoints — inputs: none — outputs: none — `bin`
- `generate_kinases_from_biomart.R` — Kinases wrapper that runs the generic BioMart generator (InterPro + GO), writes temp outputs and union — inputs: Manning file, BioMart — outputs: `genesets/curated/kinases/outputs/human_union.csv` and temp files (placed in `output/temp/` as `human_domain_interpro.csv` and `human_go_filtered.csv`) — `bin`
- `kinase_validation.R` — Alternate/duplicate validation script (similar to comprehensive_kinase_validation.R) — inputs: curated kinases CSV — outputs: `mouse_kinome_validation_results.csv` — `bin` or `archive` (duplicate)
- `map_human_to_mouse_uniprot.R` — Map human kinases to mouse via UniProt/BioMart — inputs: human kinases lists — outputs: mapping CSVs — `lib`
- `map_mouse_uniprot_biomart.R` — Map mouse UniProt to BioMart IDs — inputs: UniProt mapping — outputs: mapped CSVs — `lib`
- `mappings_index.R` — Index of mapping helper scripts — inputs: none — outputs: none — `bin`
- `merge_kinase_uniprot_validation.R` — Merge UniProt validation into kinases baseline — inputs: uniprot mapping files — outputs: merged CSV — `lib`
- `scan_invalid_entrez.R` — Scan kinases list for invalid Entrez IDs — inputs: kinases CSV — outputs: report CSV — `lib`
- `source_shiny_app.R` — Entry to run Shiny app (if present) — inputs: none — outputs: runs app — `bin`
- `summarize_kinase_groups.R` — Summarize Manning/HGNC groups and produce tables — inputs: kinases CSV — outputs: summary tables/CSV — `lib`
- `validate_against_genome.R` — Validate gene presence against genome reference files — inputs: genome files, kinases CSV — outputs: validation report — `lib`
- `validations_index.R` — Index to run validation scripts — inputs: none — outputs: none — `bin`

Recommendations

- Create directories: `scripts/kinases/bin`, `scripts/kinases/lib`, `scripts/kinases/utils`, `scripts/kinases/archive` and keep `scripts/kinases/temp` for intermediates.
- Move/organize scripts as recommended above. Keep duplicate validation scripts in `archive/` until canonicalized.
- Add `scripts/kinases/README.md` and `scripts/kinases/kinases_pipeline.R` as a lightweight runner to execute fetch -> generate -> validate steps.

Inventory file written by the AI agent on: 2025-12-28

Note: On 2025-12-28 the original top-level scripts were moved into `scripts/kinases/archive/` to make `bin/` and `lib/` the canonical locations. See `scripts/kinases/archive/` for preserved originals.

Checksums & temp policy (updated):

- MD5 files are colocated with their parent CSVs (no separate `checksums/` folder). Example: `genesets/curated/kinases/kinases_human_union.csv.md5`.
- Intermediate/temp CSVs are now stored in `output/temp/` by default; their `.md5` checksum files are kept in the same folder (example: `output/temp/human_domain_interpro.csv.md5`).

Archived files (moved):
- annotate_lipid_kinases_from_kegg.R
- annotate_metabolic_kinases_from_kegg.R
- augment_matching_with_aliases.R
- build_kinome_annotation.R
- clean_comprehensive_mouse_kinome.R
- comprehensive_kinase_validation.R
- export_kinase_gmt.R
- fetch_comprehensive_mouse_kinome_biomart.R
- fetch_definitive_mouse_kinome.R
- fetch_hgnc_kinase_groups.R
- fetch_kinhub_and_merge.R
- fetch_kinome.R
- fetch_mouse_kinome_from_kinhub.R
- fetch_uniprot_kinase_mapping_api.R
- fetch_uniprot_mouse_mapping_api.R
- fetch_val_sources_and_merge.R
- fetchers_index.R
- generate_kinases_from_biomart.R
- kinase_validation.R
- kinases_pipeline.R
- map_human_to_mouse_uniprot.R
- map_mouse_uniprot_biomart.R
- mappings_index.R
- merge_kinase_uniprot_validation.R
- scan_invalid_entrez.R
- source_shiny_app.R
- summarize_kinase_groups.R
- validate_against_genome.R
- validations_index.R
