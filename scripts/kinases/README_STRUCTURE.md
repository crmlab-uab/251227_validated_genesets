# Canonical kinases pipeline directory structure (post-cleanup)

- scripts/kinases/
    - 01_fetch_geneset_BioMart.R
    - 02_fetch_validation_sources.R
    - 03_build_annotations.R
    - 04_map_human_to_mouse.R
    - 05_export_gmt.R
    - run_kinases_pipeline.sh
    - annotations_index.R
    - PIPELINE.md
    - README.md
    - INVENTORY.md
    - archives/ (legacy or backup scripts)
    - lib/ (all helper/library R scripts)

- All outputs: curated/kinases/outputs/
- All inputs: curated/kinases/inputs/

## Rules
- No scripts or outputs in bin/ or nested folders
- All IO is config-driven (YAML)
- No .DS_Store or duplicate logs

---

**This structure is now enforced. All future scripts and outputs must follow this layout.**
