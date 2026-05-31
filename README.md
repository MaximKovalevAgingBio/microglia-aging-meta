# Cross-Species and Cross-Condition Meta-Analysis of Microglial Dynamics

This repository contains the official computational pipeline and source code for the single-nucleus meta-analysis of microglial alterations in brain aging and neurodegenerative diseases (Alzheimer's disease, Parkinson's disease, amyotrophic lateral sclerosis, frontotemporal lobar degeneration) across humans and non-human primates (NHPs).

## Data Availability
Our analysis is built entirely on open-access data from GEO and other resources, as described in the original articles. The pre-processed versions of our datasets, including information about donors, predicted cell types, and clusters, are available on Zenodo at https://doi.org/10.5281/zenodo.20445437.

## Repository Structure
* `01_dataset1_processing.ipynb` — Reference preprocessing, quality control filtration, and Leiden clustering pipeline, followed by marker gene selection, quasi-binomial logistic regression, and pseudobulk analysis (demonstrated on Dataset 1). Optimized for Google Colab execution. [HTML](https://htmlpreview.github.io/?https://github.com/MaximKovalevAgingBio/microglia-aging-meta/blob/main/01_dataset1_processing.html)
* `02_human_vs_nhp_pipeline.ipynb` — Comparative evolutionary pipeline running Linear Mixed-Effects Model (LMM) regression and donor-specific logFC analysis between human and macaque datasets. [HTML](https://htmlpreview.github.io/?https://github.com/MaximKovalevAgingBio/microglia-aging-meta/blob/main/02_human_vs_nhp_pipeline.html)
* `03_cscore_coexpression_matrices.ipynb` — Contains analytical code to evaluate gene co-expression values for 183 pre-selected genes based on the CS-CORE method. [HTML](https://htmlpreview.github.io/?https://github.com/MaximKovalevAgingBio/microglia-aging-meta/blob/main/03_cscore_coexpression_matrices.html)
* `04_pseudobulk_analysis.ipynb` — Contains analytical code to evaluate cross-dataset pseudobulk gene expression behavior across aging and disease conditions. [HTML](https://htmlpreview.github.io/?https://github.com/MaximKovalevAgingBio/microglia-aging-meta/blob/main/04_pseudobulk_analysis.html)
* `05_tf_target_networks.ipynb` — Analytical framework for processing the pan-dataset TF-target interaction matrix derived from hdWGCNA to evaluate target gene regulation profiles. [HTML](https://htmlpreview.github.io/?https://github.com/MaximKovalevAgingBio/microglia-aging-meta/blob/main/05_tf_target_networks.html)
* `06_epigenetics_linger.ipynb` — Chromatin accessibility profiling (snATAC-seq) and cis-regulatory element inference via LINGER.
* `R1_hdWGCNA_commands.R` — Executable pipeline for hdWGCNA execution (modules construction, followed by TF-target network creation).
* `R2_CSCORE_commands.R` — CS-CORE pipeline to calculate all-against-all gene co-expression matrix in the Seurat object. Example for Dataset 1 given.

## Pipeline Generalization Note
For `01_dataset1_processing.ipynb` and `R2_CSCORE_commands.R`.
The remaining integrated datasets were processed using the exact same logical pipelines as demonstrated in the reference files, with minor adjustments applied for dataset-specific metadata column names (e.g., sample vs. sample_region) or model configurations (e.g., Bonferroni correction for DS3, species-specific model for macaque datasets).

## Citation
*Manuscript in preparation*
