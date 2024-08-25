# Breast Cancer Prognostic Analysis

This repository contains R scripts and supporting files for analyzing breast cancer prognostic signatures, with a particular focus on older patients (aged 70 years and older). The primary goal is to assess the prognostic capacity of various gene signatures to aid in treatment decision-making.

## Project Structure

- **main_signatures.R**: This script contains the implementation of the gene signature analysis. It utilizes several R packages for data manipulation, statistical analysis, and visualization.
- **survival_analysis.R**: This script focuses on survival analysis, applying Kaplan-Meier and Cox proportional hazards models to evaluate the prognostic value of gene signatures.

## Benchmarking Gene Signatures

This project benchmarks the prognostic capacity of seven different gene signatures across 39 breast cancer datasets, consisting of 9,583 patients. After filtering for patients aged 70 years or older, with estrogen receptor (ER) positivity and available survival data, the analysis was conducted on 871 patients. The signatures tested include:

1. **Genomic Grade Index (GGI)**
2. **70-Gene Signature**
3. **Recurrence Score (RS)**
4. **Cell Cycle Score (CCS)**
5. **PAM50 Risk-of-Recurrence Proliferation (ROR-P)**
6. **PAM50 Signature**
7. **[Your Seventh Signature]** (If applicable, add here)

### Study Subgroups
- **All Patients**: All 871 patients were analyzed.
- **ER-Positive/Lymph Node-Positive (ER+/LN+)**: 335 patients.
- **ER-Positive/Lymph Node-Negative (ER+/LN-)**: 374 patients.

### Analysis Methods
- **Kaplan-Meier Analysis**: Survival curves were generated to compare the prognostic capacity of each signature.
- **Cox Proportional Hazards Modeling**: Multivariable models were used to assess the independent prognostic value of each signature.

## Setting Up the Environment

To ensure a consistent and reproducible environment for running the R scripts, a conda environment file (`environment.yml`) is provided. This file includes all the necessary packages and their specific versions, making it easy to set up the environment on any system.

### Steps to Create the Conda Environment

1. **Clone the Repository:**

   Begin by cloning this repository to your local machine:

   ```bash
   git clone https://github.com/yourusername/your-repo-name.git
   cd your-repo-name
    ```
2. **Execute the following to generate a conda environment with all you need to execute the scripts:**
   ```bash
    conda env create -f environment_benchmark.yml
    ```
2. **Once you have created a new conda environment, make sure to activate it before executing the scripts:**
   ```bash
    conda activate elderly_env
    ```
