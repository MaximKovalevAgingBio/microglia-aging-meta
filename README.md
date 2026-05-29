# Cross-Species and Cross-Condition Meta-Analysis of Microglial Dynamics

This repository contains the official computational pipeline and source code for the single-nucleus meta-analysis of microglial alterations in brain aging and neurodegenerative diseases (Alzheimer's disease, Parkinson's disease, amyotrophic lateral sclerosis, frontotemporal lobar degeneration) across humans and non-human primates (NHPs).

## Data Availability
Our analysis is built entirely on open-access data from GEO and other resources, as described in the original articles. The pre-processed versions of our datasets, including information about donors, predicted cell types, and clusters, are available on Zenodo at https://doi.org/10.5281/zenodo.20445437.

## Repository Structure
* `01_dataset1_processing.ipynb` — Reference preprocessing, quality control filtration, and Leiden clustering pipeline (demonstrated on Dataset 1). Optimized for Google Colab execution.
* `02_Human_vs_NHP_pipeline.ipynb` — Comparative evolutionary pipeline running Linear Mixed-Effects Model (LMM) regression and donor-specific logFC analysis between human and macaque datasets.
* `03_downstream_analysis.ipynb` — Downstream pipeline integration for CS-CORE co-expression outputs, transcription factor (TF)-target networks, and pseudo-bulk processing.
* `04_epigenetics_linger.ipynb` — Chromatin accessibility profiling (snATAC-seq) and cis-regulatory element inference via LINGER.
* `hdWGCNA_commands.txt` — Terminal guide containing console commands and parameter layouts used to execute hdWGCNA modules.
* `CS-CORE_commands.txt` — Terminal guide containing console commands for CS-CORE matrix calculation.

## Pipeline Generalization Note
The remaining integrated datasets were processed using the exact same logical pipelines as demonstrated in the reference files, with minor adjustments applied for dataset-specific metadata column names (e.g., sample vs. sample_region) or model configurations (e.g., Bonferroni correction for DS3, disease-specific models for macaque tissue).

## Citation
*Manuscript in preparation*
