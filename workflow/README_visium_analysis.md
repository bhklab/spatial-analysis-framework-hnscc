# Code overview (Parts 1–4)

This folder contains the analysis workflow scripts organized by manuscript stages (Part 1–4).
The repository is intended for methodological transparency; running end-to-end requires access to the corresponding input objects/files described in each Part README.

## Folder structure

visium_analysis/
├── 01_part1/ # Preprocessing + spot/cell annotations + deconvolution + malignant calling
├── 02_part2/ # Malignant/non-malignant characterization + neighborhood analysis
├── 03_part3/ # inferCNV + CellChat (compute-focused; minimal/no visualization)
└── 04_part4/ # Figure generation (visualization only)

## 01_part1 — Preprocessing & core annotations

Goal: Build and annotate spatial Seurat objects (pathology/QC/deconvolution/malignant calling).

Scripts:
- 01_create_seurat_objects.R — create Seurat objects from raw/processed inputs
- 02_add_pathologist_annotations.R — add pathology annotations to spots
- 03_qc_and_spotsweeper.R — QC filtering + SpotSweeper
- 04_build_sc_reference.R — build single-cell reference object(s)
- 05_deconvolution.R — spatial deconvolution (e.g., RCTD/CARD as used)
- 06_export_for_scmalignantfinder.R — export inputs for scMalignantFinder
- 07_scMalignant.ipynb — scMalignantFinder execution (notebook)
- 08_import_scmalignant_results.R — import malignancy predictions back into Seurat
- 09_malignant_annotation.R — finalize malignant vs non-malignant labels

Documentation:
- README_part1.md

## 02_part2 — Characterization & neighborhood analyses

Goal: Characterize malignant/non-malignant programs and spatial neighborhoods.

Scripts:
- 01a_global_integration_and_umap.R — global integration + UMAP/cluster framework
- 01b_malignant_signatures_and_subclustering.R — malignant signatures + subclustering framework (labels refined later in Part 4)
- 02_nonmalignant_spot_characterization.R — non-malignant state characterization
- 03_neighborhood_analysis.R — neighborhood feature extraction/analysis
- 04_identifying_near_far.R — define near vs far regions relative to malignant interface
- 05_dge_near_vs_far.R — differential expression (near vs far)

Documentation:
- README_part2.md

## 03_part3 — inferCNV + CellChat (compute)

Goal: inferCNV and cell–cell communication analyses (compute outputs for later visualization).

Scripts:
- 01_infercnv_round1_round2.R — inferCNV Round 1 (malignant vs non-malignant) + Round 2 (malignant subclusters)
- 02_cellchat.R — CellChat analysis (Tier1/Tier2/Tier2b/Tier2c)

Documentation:
- README_part3.md

## 04_part4 — Figure generation (visualization only)

Goal: Generate manuscript and supplementary figures from Part 1–3 outputs.

Files:
- 00_common.R — shared plotting helpers/theme/palette conventions
- fig1.R … fig5.R — main figures
- suppl_fig1.R, suppl_fig3.R, suppl_fig4.R — supplementary figures

Documentation:
- README_part4.md

## Software environment (packages + versions)

The following R packages (and versions) were used in this analysis:

- Seurat (v.5.3.1)
- DropletUtils (v.1.28.1)
- yaml (v.2.3.11)
- RCTD (v.2.2.0)
- SpotSweeper (v.1.4.0)
- SpatialExperiment (v.1.18.1)
- SummarizedExperiment (v.1.38.1)
- scuttle (v.1.18.0)
- dplyr (v.1.1.4)
- Matrix (v.1.7-4)
- CARD (v.1.0.0)
- stringr (v.1.6.0)
- future (v.1.68.0)
- AUCell (v.1.32.0)
- readxl (v.1.4.5)
- clustree (v.0.5.1)
- purr (v.1.2.0)
- readr (v.2.1.6)
- spdep (v.1.4-1)
- tidyr (v.1.3.1)
- igraph (v.2.2.1)
- tibble (v.3.3.0)
- infercnv (v.1.22.0)
- CellChat (v.1.6.1)
- patchwork (v.1.3.2)
- ggplot2 (v.4.0.1)
- BiocGenerics (v.0.54.1)
- BiocManager (v.1.30.27)
- Nebulosa (v.1.18.0)
- RColorBrewer (v.1.1-3)
- SCpubr (v.3.0.0)
- ape (v.5.8-1)
- cowplot (v.1.2.0)
- ggtree (v.3.16.3)
- grid (v.4.5.1)
- phangorn (v.2.12.1)
- pheatmap (v.1.0.13)
- presto (v.1.0.13)
- treeio (v.1.32.0)
- GSVA (v.2.2.1)
- circlize (v.0.4.17)
- forcats (v.1.0.1)
- clusterProfiler (v.4.16.0)
- msigdbr (v.25.1.1)
- org.Hs.eg.db (v.3.12.0)
- writexl (v.1.5.4)

## Notes

- Large intermediate objects (e.g., *.rds) and generated figures are typically not committed to GitHub.
- Each Part README documents expected inputs/outputs and key steps.
