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

### Study Subgroups
- **All Patients**: All 871 patients were analyzed.
- **ER-Positive/Lymph Node-Positive (ER+/LN+)**: 335 patients.
- **ER-Positive/Lymph Node-Negative (ER+/LN-)**: 374 patients.

### Analysis Methods
- **Kaplan-Meier Analysis**: Survival curves were generated to compare the prognostic capacity of each signature.
- **Cox Proportional Hazards Modeling**: Multivariable models were used to assess the independent prognostic value of each signature.


### To extract the data please go to this [bitbucket repository](https://bitbucket.org/tobingroup/elderly/src/master/)

## Setting Up the Environment

To ensure a consistent and reproducible environment for running the R scripts, a conda environment file (`environment.yml`) is provided. This file includes all the necessary packages and their specific versions, making it easy to set up the environment on any system.

### Steps to Create the Conda Environment

1. **Clone the Repository:**

   Begin by cloning this repository to your local machine:

   ```bash
   git clone https://github.com/MiguelCastresana/elderly.git
   cd elderly
    ```
2. **Execute the following to generate a conda environment with all you need to execute the scripts:**
   ```bash
    conda env create -f environment.yml
    ```
3. **Once you have created a new conda environment, make sure to activate it before executing the scripts:**
   ```bash
    conda activate elderly_env
    ```
4. **Download the data necessary to run these study, [here](https://drive.google.com/drive/folders/1KkRhLCEQdkR4TjqPWwB2A9-akbrVgyrF?usp=sharing). Important to add the downloaded folder inside the cloned repository**

### Run the Main Gene Signature Analysis

The main gene signature analysis is performed by running the main_signatures.R script. This script processes the data, applies the gene signatures, and evaluates their prognostic capacity across different patient subgroups.

To run the analysis, execute the following command in your terminal:
   ```bash
    Rscript main_signatures.R
   ```
### Run the Survival Analysis

The survival analysis is performed using the survival_analysis.R script. This script generates Kaplan-Meier survival curves and performs multivariable Cox proportional hazards modeling to assess the independent prognostic value of each gene signature.

To run the survival analysis, execute the following command in your terminal:

   ```bash
    Rscript survival_analysis.R
   ```

### Results

Running these scripts will produce the following outputs:

1. **Kaplan-Meier Survival Curves:**
Visual Outputs: The survival_analysis.R script generates Kaplan-Meier survival curves for each gene signature across different patient subgroups (all patients, ER+/LN+, ER+/LN-). These plots will be saved as image files (e.g., .png or .pdf) in the specified output directory.
Interpretation: These curves show the probability of survival over time for different risk groups as classified by each gene signature.

2. **Cox Proportional Hazards Models:**
Tabular Outputs: The script will output tables summarizing the results of the multivariable Cox proportional hazards models. These tables include hazard ratios (HR), confidence intervals (CI), and p-values for each gene signature within the different patient subgroups.
Interpretation: These results help to understand the independent prognostic value of each gene signature when other variables (e.g., age, tumor grade) are taken into account.

3. **Gene Signature Scores:**
Data Outputs: The main_signatures.R script outputs the calculated scores for each gene signature for every patient in the dataset. These scores are saved in data files (e.g., .csv format) for further analysis or validation.
Interpretation: These scores indicate the risk level predicted by each gene signature for individual patients, which can be correlated with survival outcomes.


**The publication can be found, [here](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-024-01797-7) **


**Contact**
Miguel Castresana Aguirre (miguel.castresana.aguirre@ki.se)
