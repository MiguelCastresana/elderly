# Legacy Data-Preparation Scripts

The main reproducible workflow lives in the top-level `R/` scripts.

The `R/legacy/` folder contains older data-ingestion and mapping scripts that were useful during study assembly but still reference the original local raw-data export layout. They are kept for provenance, but they are not required for the main `R/main_signatures.R` workflow when the shared `data_bitbucket/` folder is already available.

Legacy scripts:

- `R/legacy/datasets_cellcycle.R`
- `R/legacy/fetch_parse_annotation_GEX_data.R`
- `R/legacy/mapping.R`

