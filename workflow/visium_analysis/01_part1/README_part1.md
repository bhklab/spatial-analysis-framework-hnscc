# Part 1 — Spatial preprocessing, reference building, deconvolution, and malignant calling

This folder contains the end-to-end preprocessing pipeline used to generate analysis-ready Visium Seurat objects, build single-cell references, run deconvolution, and assign malignant vs non-malignant spots (including scMalignantFinder).

## Overview and run order

Run scripts in this order:

1. `01_create_seurat_objects.R`  
2. `02_add_pathologist_annotations.R`  
3. `03_qc_and_spotsweeper.R`  
4. `04_build_sc_reference.R`  
5. `05_deconvolution.R`  
6. `06_export_for_scmalignantfinder.R`  
7. `07_run_scmalignantfinder.ipynb`  
8. `08_import_scmalignant_results.R`  
9. `09_malignant_annotation.R`

## Inputs (what you need before running)

### A) SpaceRanger outputs (Visium)
For each spatial sample you need the SpaceRanger output folder containing:
- `filtered_feature_bc_matrix/` **or** `filtered_feature_bc_matrix.h5`
- `spatial/` (tissue image + coordinates, etc.)

### B) Sample sheet(s)
You should have a sample manifest (CSV) that maps each `sample_id` to the corresponding SpaceRanger folder name.

Recommended columns (minimal):
- `sample_id`
- `spaceranger_dir`

Optional but useful:
- `timepoint` (Pre/Post)
- `patient_id`

### C) Pathologist annotation CSV(s) (per sample)
For each sample, a CSV file named:
- `<sample_id>.csv`

Required columns:
- `Barcode`
- `Patho_Anno_Jinsu` (or update script setting if your column name differs)

### D) Single-cell reference data (for deconvolution)
This pipeline supports building references from multiple sources (e.g., Choi, Quah, Puram).
You will need the raw reference input files used by `04_build_sc_reference.R` (paths are defined in that script).

### E) scMalignantFinder
You need a working scMalignantFinder installation/environment for running the notebook:
- `07_run_scmalignantfinder.ipynb`

You will also need the scMalignantFinder pretrained model directory (pretrain_dir) if required by your setup.

### F) Clinical metadata for final annotation
For `09_malignant_annotation.R` you need:
- `sample_info.csv`
- `patient_info.csv`

These are merged onto the Seurat metadata.

## Outputs (what each step produces)

All paths below are conceptual; adjust to your project’s `output/` location.

### Step 1 — Build Seurat objects
**Script:** `01_create_seurat_objects.R`  
**Input:** SpaceRanger output folders  
**Output:** `output/part1/01_seurat_raw/<sample_id>.rds`  
**Notes:** Ensures `filtered_feature_bc_matrix.h5` exists (creates from `filtered_feature_bc_matrix/` when needed).

### Step 2 — Add pathologist annotation
**Script:** `02_add_pathologist_annotations.R`  
**Input:** `output/part1/01_seurat_raw/*.rds` + per-sample annotation CSV  
**Output:** `output/part1/02_pathologist_annot/<sample_id>.rds`  
**Also writes:** `pathologist_annotation_summary.csv` (barcode match summary)

### Step 3 — QC + SpotSweeper filtering
**Script:** `03_qc_and_spotsweeper.R`  
**Input:** `output/part1/02_pathologist_annot/*.rds`  
**Output:** `output/part1/03_qc_spotsweeper/<sample_id>.rds`  
**Also writes:**
- `SpotSweeper_outlier_summary.csv`
- `Supplementary_SpotFiltering_Summary.csv`

Filtering logic:
- remove local outliers flagged by SpotSweeper (`sum`, `detected`, `mito%`)
- apply global thresholds (default in script): `nCount_Spatial > 200` and `nFeature_Spatial > 200`

### Step 4 — Build single-cell references
**Script:** `04_build_sc_reference.R`  
**Input:** reference datasets (Choi / Quah / Puram; see script for file paths)  
**Output:** `output/part1/04_sc_reference/*.rds`  
**Also writes:** diet (memory-light) versions under `output/part1/04_sc_reference/refs_diet/`

### Step 5 — Deconvolution (CARD)
**Script:** `05_deconvolution.R`  
**Input:** filtered Visium objects + diet references  
**Output:** `output/part1/04_deconvolution/<ref>/<sample_id>_CARD.rds`  
**Also writes:** per-sample proportion tables if available:
- `<sample_id>_CARD_proportions.csv`

### Step 6 — Export data for scMalignantFinder
**Script:** `06_export_for_scmalignantfinder.R`  
**Input:** deconvolved (or latest) per-sample Seurat objects  
**Output:** `output/part1/05_scmf_input/`
- `<sample_id>_counts.csv` (genes x spots)
- `<sample_id>_metadata.csv` (spots x metadata)

### Step 7 — Run scMalignantFinder (Python)
**File:** `07_run_scmalignantfinder.ipynb`  
**Input:** `output/part1/05_scmf_input/`  
**Output:** `output/part1/06_scmf_output/`
- `<sample_id>_malignant_results.csv` (spot-level predictions)
- optionally: `*_with_malignancy.h5ad` (if saved)

### Step 8 — Import scMalignantFinder results into Seurat
**Script:** `08_import_scmalignant_results.R`  
**Input:** Seurat objects + `output/part1/06_scmf_output/*_malignant_results.csv`  
**Output:** `output/part1/07_malignant_annot/<sample_id>_Annotated.rds`

Adds metadata columns:
- `malignancy_probability`
- `scMalignantFinder_prediction`

### Step 9 — Final malignant call + broad cell type label
**Script:** `09_malignant_annotation.R`  
**Input:** annotated objects + clinical metadata CSVs  
**Output:**
- `output/part1/08_final_malignant_annotation/AllST.rds`
- `output/part1/08_final_malignant_annotation/per_sample/<sample_id>.rds`

Final malignant rule (as implemented in the script):
- Pathologist annotation indicates SCC
- scMalignantFinder predicts malignant
- deconvolution consensus (dominant malignant calls across methods) meets threshold

Also derives a broad non-malignant label (`Broad_CT`) from deconvolution proportions.

## Reproducibility notes

- Avoid committing large objects (RDS, figures) to GitHub.
- Keep personal machine paths out of committed code when possible.
- If using a local `config.yml`, commit only a `config.example.yml` template and ignore the real config in `.gitignore`.

## Troubleshooting

- If barcode matching in Step 2 yields near-zero matches, confirm:
  - annotation CSV barcode format (`AAAC...` vs `AAAC...-1`)
  - Seurat barcode format (`<sample>_AAAC...-1` vs `AAAC...-1`)
- If SpotSweeper flags most spots, check:
  - mitochondrial gene prefix (`^MT-`)
  - whether counts slot/assay is correct (`Spatial`)
- If scMalignantFinder results do not import, confirm:
  - result CSV rownames are spot barcodes matching Seurat `colnames()`

