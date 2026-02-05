# Raw Data Directory

## Purpose

This directory is reserved for **immutable raw data** files used in the spatial transcriptomics analyses presented in this study.  
Raw data are **not tracked in this repository** due to patient privacy, consent restrictions, and data-sharing agreements.

All analyses in this repository assume that raw data have been obtained separately and placed in user-defined local paths.

---
## Data Access Instructions

**No raw data files are included in this repository.**  
Access to raw data depends on the dataset type and is governed by ethical approvals and data-use agreements.

This study integrates private spatial transcriptomics data and public/private bulk RNA-seq datasets, along with curated immuno-oncology gene expression signatures, for external validation.

---

## Spatial Transcriptomics Data (Discovery Cohort)

### INSPIRE and IO-KIN Clinical Trials

Spatial transcriptomic profiling was performed on tumor biopsies from immune checkpoint blockade (ICB)–naïve patients with recurrent or metastatic head and neck squamous cell carcinoma (RM-HNSCC) enrolled in two investigator-initiated phase II trials:

- **INSPIRE** [NCT02644369](https://clinicaltrials.gov/study/NCT02644369)  
- **IO-KIN** [NCT04606940](https://clinicaltrials.gov/study/NCT04606940)

**Data type**
- 10x Genomics Visium spatial transcriptomics from FFPE tissue sections (10 μm)  
- Cell Ranger–generated outputs including gene expression matrices, spatial coordinates, and histology images  

**Access**
- Data are not publicly available due to patient privacy and institutional restrictions  
- Access requires institutional approval and appropriate data-sharing agreements; interested users should contact the corresponding authors  

---

## Single-Cell RNA-seq Reference Datasets

Public single-cell RNA-seq datasets were used exclusively as reference atlases for spatial deconvolution and cell-type annotation.

**Referenced datasets**
- [GSE181919](https://pubmed.ncbi.nlm.nih.gov/36828832/)  
- [GSE182227](https://pubmed.ncbi.nlm.nih.gov/37012457/)  
- [GSE188737](https://pubmed.ncbi.nlm.nih.gov/36973261/)  

**Data type**
- Whole-transcriptome single-cell RNA-seq, publicly available via GEO  

These datasets are not redistributed in this repository and must be downloaded directly from their original sources.

---

## Bulk RNA-seq Data (External Validation Cohorts)

To evaluate the generalizability of spatially informed signatures, bulk RNA-seq datasets from multiple immune checkpoint blockade (ICB) cohorts were analyzed.

### ORCESTRA Platform

Bulk RNA-seq datasets used for external validation were accessed via [**ORCESTRA**](https://www.orcestra.ca), a reproducible biomedical data platform.

Please download the dataset using the following curated release:

➡️ **https://www.orcestra.ca/clinical_icb**

This dataset includes:
- Pre-treatment normalized RNA-seq expression data (TPM)  
- Clinical annotations for immune checkpoint blockade (ICB) studies  
- Associated metadata and processing documentation  

Only a subset of cohorts from ORCESTRA was included, based on predefined inclusion criteria.

---

## Gene Signature Resources

Gene expression signatures used for scoring, enrichment, and benchmarking were obtained from publicly available resources:

- **IO Signatures** — [bhklab/SignatureSets](https://github.com/bhklab/SignatureSets)

Each signature is:
- Annotated with source publication  
- Categorized (e.g., ICB-sensitive, ICB-resistant)  
- Used for GSVA, weighted mean, ssGSEA, or method-specific scoring  

---

## Inclusion & Exclusion Criteria (Bulk Validation Cohorts)

Raw datasets and gene signatures were selected based on:
- Availability of pre-treatment RNA-seq data  
- Available clinical outcomes (PFS and/or objective response)  
- PD-1 or PD-L1 monotherapy  
- Sufficient sample size per cohort  

Datasets were excluded if they were generated using targeted gene panels, involved combination therapy regimens, consisted only of post-treatment samples, or lacked sufficient clinical annotation.

Refer to the **Materials and Methods** section of the manuscript for detailed criteria.

---

## Additional Notes

- No preprocessing or normalization is performed directly on raw data within this repository.  
- All downstream processing occurs in analysis-specific scripts and intermediate data directories.  
- Users are responsible for complying with all original data-use agreements and ethical approvals.  
- This directory is **read-only** during pipeline execution; all transformations occur downstream in `data/procdata/`.  
