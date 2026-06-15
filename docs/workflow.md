# Workflow

This repository contains the analysis code for benchmarking breast cancer prognostic signatures in older patients.

## Main Entry Points

- `R/main_signatures.R`: orchestrates signature scoring and downstream analyses.
- `R/survival_analysis.R`: generates Kaplan-Meier plots and Cox model summaries after `results_final` has been prepared.
- `R/final_analysis_all_above70.R`: merges and filters results for the main >=70 cohort.
- `R/final_analysis_all_55_65.R`: runs the comparator 55-65 cohort analysis.

## Supporting Signature Scripts

- `R/oncotype_function.R`
- `R/mammaprint_function.R`
- `R/ggi_function.R`
- `R/cell_cycle_function.R`
- `R/PAM50_RORP_function.R`
- `R/PAM50_RORP_function_nomontecarlo.R`

## Legacy Scripts

Older raw-data preparation scripts are kept in `R/legacy/` for provenance. They are not required for the main analysis if `data_bitbucket/` has already been downloaded. See [legacy-scripts.md](legacy-scripts.md).

## Suggested Run Order

From the repository root:

```bash
conda env create -f environment.yml
conda activate elderly_env
Rscript R/main_signatures.R
```

The full workflow requires `data_bitbucket/` to be present.
