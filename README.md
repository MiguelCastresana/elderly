# Breast Cancer Prognostic Signatures in Older Patients

This repository contains the R analysis code for benchmarking breast cancer prognostic signatures in older patients, with a focus on patients aged 70 years and older.

The study evaluates whether established gene-expression signatures retain prognostic value in older, estrogen receptor-positive breast cancer cohorts and whether they can support treatment decision-making.

Publication: [Benchmarking breast cancer prognostic signatures for elderly patients](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-024-01797-7)

## Analysis Scope

The analysis benchmarks prognostic signatures across public breast cancer cohorts and evaluates survival associations in clinically relevant patient subgroups.

Signatures included:

- Genomic Grade Index (GGI)
- 70-gene signature / MammaPrint
- Recurrence Score / Oncotype DX-like score
- Cell Cycle Score (CCS)
- PAM50 subtype
- PAM50 Risk of Recurrence with proliferation (ROR-P)

Main cohorts:

- All eligible older patients
- ER-positive / lymph node-positive patients
- ER-positive / lymph node-negative patients
- Comparator 55-65 age group

## Repository Layout

```text
.
├── R/                  # Analysis scripts and signature implementations
├── docs/               # Workflow and data notes
├── tools/              # Lightweight project checks
├── environment.yml     # Conda environment
└── README.md
```

Important entry points:

- `R/main_signatures.R`: main workflow for signature scoring, result merging, and survival analysis.
- `R/survival_analysis.R`: Kaplan-Meier and Cox proportional hazards analyses. This expects merged results prepared by the main workflow.
- `R/final_analysis_all_above70.R`: main older-patient cohort merge and filtering.
- `R/final_analysis_all_55_65.R`: comparator cohort merge and filtering.

See [docs/workflow.md](docs/workflow.md) for more detail.

Older raw-data preparation scripts are kept in `R/legacy/` for provenance. They are not part of the main run once `data_bitbucket/` has been downloaded.

## Data

The analysis requires external study data that is not committed to this repository.

Download the data folder from the shared link:

<https://drive.google.com/drive/folders/1KkRhLCEQdkR4TjqPWwB2A9-akbrVgyrF?usp=sharing>

Place the downloaded folder in the repository root and name it:

```text
data_bitbucket/
```

See [docs/data.md](docs/data.md) for the expected file structure.

## Environment

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate elderly_env
```

## Quick Check

Before running the full analysis, verify that the repository structure and R syntax are valid:

```bash
Rscript tools/check-project.R
```

This check does not require the private data folder. If `data_bitbucket/` is missing, it reports that the full analysis is skipped.

## Run The Analysis

From the repository root:

```bash
Rscript R/main_signatures.R
```

The workflow reads inputs from `data_bitbucket/` and writes generated result files and plots under `data_bitbucket/final_results/`.

## Outputs

The scripts generate:

- Signature scores for each eligible sample.
- Merged analysis tables for older and comparator cohorts.
- Kaplan-Meier survival plots.
- Cox proportional hazards model summaries.
- Subgroup analyses by ER and lymph-node status.

## Notes

- The repository now uses portable project-root detection through `R/paths.R`.
- The code assumes the downloaded data folder keeps the original `data_bitbucket/` name.
- Generated data, plots, and local R session files are ignored by Git.

## Contact

Miguel Castresana Aguirre  
[miguel.castresana.aguirre@ki.se](mailto:miguel.castresana.aguirre@ki.se)
