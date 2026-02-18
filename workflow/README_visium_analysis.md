# Scripts Directory — Visium Spatial Analysis

## Purpose

This directory provides modular R scripts implementing a comprehensive spatial transcriptomics analysis workflow for 10x Visium datasets.

The workflow includes:

- Spatial object construction and QC
- Pathology-guided spot annotation
- Deconvolution and malignant calling
- Malignant and non-malignant state characterization
- Spatial neighborhood analysis
- inferCNV and cell–cell communication modeling
- Manuscript-quality visualization

The pipeline is modular and structured for reproducibility. 

---

## Key Capabilities

- Construction and annotation of spatial Seurat objects
- Pathology-integrated malignant classification
- Spatial deconvolution (RCTD / CARD)
- Neighborhood feature extraction and interface analysis
- inferCNV-based genomic instability analysis
- CellChat-based ligand–receptor communication analysis
- Differential expression (near vs far interface)
- Automated manuscript figure generation

---

# Modular Scripts (`workflow/visium_analysis/`)

The workflow is organized into eight functional modules aligned with the repository-wide framework.

---

## (1) runProcData

Preprocessing and core spatial annotations.

Includes:
- `01_create_seurat_objects.R`
- `02_add_pathologist_annotations.R`
- `03_qc_and_spotsweeper.R`
- `04_build_sc_reference.R`
- `05_deconvolution.R`
- `06_export_for_scmalignantfinder.R`
- `07_scMalignant.ipynb`
- `08_import_scmalignant_results.R`
- `09_malignant_annotation.R`

Outputs:
- Annotated Seurat objects
- Deconvolution matrices
- Malignant vs non-malignant labels

---

## (2) runSigScore

Computation of malignant programs and spatial signature scores.

Includes:
- Malignant signature derivation
- AUCell scoring
- Subclustering-based program definition

Primary script:
- `01b_malignant_signatures_and_subclustering.R`

---

## (3) runSigAssoc

Association and differential analyses.

Includes:
- `05_dge_near_vs_far.R`
- Malignant vs non-malignant comparisons

Outputs:
- Differential expression tables
- LogFC and adjusted p-values

---

## (4) runMeta

Cross-sample spatial integration.

Includes:
- `01a_global_integration_and_umap.R`

Performs:
- Cross-sample integration
- Dimensionality reduction
- Global clustering framework

---

## (5) runCorr

Cell–cell interaction and communication analyses.

Includes:
- `02_cellchat.R`

Outputs:
- Ligand–receptor interaction networks
- Pathway-level interaction summaries

---

## (6) runSigCluster

Clustering of malignant subpopulations and spatial states.

Includes:
- Malignant subclustering logic
- Neighborhood clustering frameworks

Associated scripts:
- `03_neighborhood_analysis.R`

---

## (7) sigDistanceFunction

Spatial proximity and interface modeling.

Includes:
- `04_identifying_near_far.R`

Computes:
- Spatial adjacency metrics
- Distance-to-interface measures
- Near vs far classification

---

## (8) runVisualization

Manuscript and supplementary figure generation.

Includes:
- `00_common.R`
- `fig1.R` – `fig5.R`
- `suppl_fig1.R`
- `suppl_fig3.R`
- `suppl_fig4.R`

Outputs:
- Publication-ready figures

---


# Input Requirements

Each Visium sample requires:

- 10x Visium output directory
- Spot coordinate files
- Histology image (H&E)
- Pathologist annotation files
- Optional single-cell reference data
- Clinical metadata (if applicable)

---

# Software Environment

The following R packages (and versions) were used:

- Seurat (v5.3.1)
- RCTD (v2.2.0)
- SpotSweeper (v1.4.0)
- CARD (v1.0.0)
- infercnv (v1.22.0)
- CellChat (v1.6.1)
- AUCell (v1.32.0)
- GSVA (v2.2.1)
- clusterProfiler (v4.16.0)
- msigdbr (v25.1.1)
- org.Hs.eg.db (v3.12.0)

Additional dependencies are documented within each script.

---

# Notes

- Large intermediate objects (.rds) and generated figures are not committed.
- Each module is self-contained and documents required inputs/outputs.
- Figure scripts assume outputs from preprocessing and analysis modules are available.
- This repository is intended for methodological transparency and reproducibility.

