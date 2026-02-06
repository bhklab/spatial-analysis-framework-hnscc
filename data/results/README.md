# Results Directory

## Purpose

This directory contains all **derived outputs** generated from the spatial transcriptomics
and external validation analyses presented in this study. These results include:

- Spatially informed gene signature scores
- Clinical association analyses (PFS and response)
- Meta-analysis summaries across cohorts
- Pathway enrichment and clustering results
- Publication-ready figures and visualizations

All outputs are generated from **preprocessed inputs** and contain **no raw sequencing data**.

---

## Directory Structure

```console
results/
├── assoc/        # Association analyses (Cox, logistic regression)
├── cluster/      # Clustering results (e.g., signature/pathway clustering)
├── cor/          # Correlation analyses between signatures or features
├── data/         # Processed cohort-level data objects used for downstream analyses
├── Fig/          # Figures (KM plots, boxplots, forest plots, volcano plots)
├── meta/         # Meta-analysis results across cohorts
├── score/        # Gene signature scores (including spatially informed signatures)
├── sig/          # Signature-specific outputs and intermediate objects
└── README.md     # Description of results structure and governance
