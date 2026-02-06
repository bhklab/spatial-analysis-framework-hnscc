# Usage Guide

This section provides step-by-step instructions to configure, run, and manage your analysis.

---

## Project Configuration

### 1. Customize Configuration Files

Edit the YAML files in the `config/` directory to reflect your study parameters:
- `config_rnaseq_io.yaml`: For processing raw input into analysis-ready format.

Each center should create a center-specific config file using the template:

```bash
config/<center_or_study_name>.yaml
# Example:
config/ICB_Gide.yaml
```

Use `config/config_local.yaml` as a reference.

---

### 2. Add Input Data

Place your raw input datasets (e.g., `.rds` files) in the `data/rawdata/` directory:

```bash
data/rawdata/ICB_Gide.rds
```

Also include gene signature metadata files:

- `signature.rda`
- `signature_information.csv`
- Precompiled `.RData` file — [Zenodo DOI: 10.5281/zenodo.18509237](https://doi.org/10.5281/zenodo.18509237)

---

## Running Your Analysis

### 1. Install Dependencies

If using [`pixi`](https://prefix.dev/docs/pixi/overview):

```bash
pixi install
```

This sets up your environment with R, Bioconductor, and required packages.

---

### 2. Run Local Processing

Prepare analysis-ready `.rda` files from raw data:

```bash
Rscript workflow/scripts/runProcData.R
```

---

### 3. Run Signature Scoring 

```bash
Rscript workflow/scripts/runSigScore.R
```

---

### 4. Run Association-analysis

```bash
Rscript workflow/scripts/runSigAssoc.R
```

### 5. Run Meta-analysis 

```bash
Rscript workflow/scripts/runMeta.R
```

---

### 6. Generate Visualizations

```bash
Rscript workflow/scripts/runVisualization.R
```

This script generates forest plots, volcano plots, and heatmaps.
**Note**: `runSigCluster.R` and `runCorr.r` are considered to study signatures' distributions. 

---

## Tips for Managing Your Data

- Use meaningful filenames (e.g., `ICB_Gide__Melanoma__PD-(L)1.rda`)