# Data Sources

## Overview

This document summarizes the **spatial transcriptomics**, **Single-cell RNA-seq**, **bulk RNA-seq**, and **RNA-based gene signature** data sources used in this study. Datasets include private clinical trial cohorts and publicly available reference and validation datasets, curated to ensure reproducibility, transparency, and consistency with the accompanying spatial transcriptomics manuscript.

---

## Data Sources

### Spatial Transcriptomics Data Sources (Discovery Cohort)

Spatial transcriptomic profiling was performed on tumor biopsies from patients enrolled in investigator-initiated immune checkpoint blockade (ICB) trials:

- **INSPIRE** [NCT02644369](https://clinicaltrials.gov/study/NCT02644369)  
- **IO-KIN** [NCT04606940](https://clinicaltrials.gov/study/NCT04606940)

**Data type**
- 10x Genomics Visium spatial transcriptomics from FFPE tissue sections (10 μm)  
- Cell Ranger–generated outputs including gene expression matrices, spatial coordinates, and histology images  

**Access**
- Restricted due to patient privacy and institutional data-use agreements  
- Data are not publicly distributed  

### Single-Cell RNA-seq Reference Data

Public single-cell RNA-seq datasets were used exclusively as reference atlases for cell-type deconvolution and spatial annotation.

- **Source**: Gene Expression Omnibus (GEO) -  [GSE181919](https://pubmed.ncbi.nlm.nih.gov/36828832/), [GSE182227](https://pubmed.ncbi.nlm.nih.gov/37012457/), [GSE188737](https://pubmed.ncbi.nlm.nih.gov/36973261/)
- **Data type**: Whole-transcriptome single-cell RNA-seq
- **Access**: Public

### Bulk RNA-seq Data (External Validation Cohorts)

#### Public Validation Datasets

Bulk RNA-seq datasets were used to evaluate the generalizability of spatially informed gene expression signatures.

- **Name**: Immune Checkpoint Blockade - RNA-Seq, and Clinical data
- **URL**: [https://www.orcestra.ca/clinical_icb](https://www.orcestra.ca/clinical_icb)
- **Access Method**: Direct download or programmatic retrieval via API (if applicable)
- **Data Format**: MultiAssayExperiment and SummarizedExperiment in R (Bioconductor)
- **Access**: Public

Only cohorts meeting predefined inclusion criteria were used.

#### Private Validation Datasets

In addition to publicly available cohorts, selected private bulk RNA-seq datasets were used for external validation of spatially informed gene expression signatures.

- **Data type**: Pre-treatment whole-transcriptome bulk RNA-seq with associated clinical outcome annotations
- **Source**: Institutional and collaborative immunotherapy cohorts not publicly released
- **Access**: Restricted under institutional review board (IRB) approvals and data-sharing agreements

These datasets were processed and analyzed using the same harmonization, normalization, and scoring procedures applied to public validation cohorts. Due to patient privacy and consent restrictions, raw data are not redistributed in this repository.

**Representative private validation cohorts include:**

| Study / Lead Author | Trial / Program            | Cancer Type                  | Data Repository | Access Status |
|--------------------|----------------------------|------------------------------|-----------------|---------------|
| Ravi et al.        | SU2C-MAR                   | Non-small cell lung cancer   | dbGaP / EGA     | Restricted    |
| McDermott et al.   | IMmotion150                | Renal cell carcinoma         | EGA             | Restricted    |
| Priestley et al.   | Hartwig Medical Foundation | Pan-cancer                   | HMF             | Restricted    |
| Rittmeyer et al.   | OAK (NCT02008227)           | Non-small cell lung cancer   | EGA             | Restricted    |
| Fehrenbacher et al.| POPLAR (NCT01903993)        | Non-small cell lung cancer   | EGA             | Restricted    |

Access to these datasets requires approval from the corresponding data custodians and compliance with institutional and ethical data-use agreements. Cohort-level inclusion, preprocessing, and harmonization procedures are described in the **Materials and Methods** and **Supplementary Data** of the manuscript.


### Gene Signature Resources

RNA-based gene expression signatures were used for spatial feature scoring, enrichment analysis, and benchmarking.

- **Name**: SignatureSets: An R Package for RNA-Based Immuno-Oncology Signatures
- **Version**: v1.0
- **URL - IO Signatures**: [bhklab/SignatureSets](https://github.com/bhklab/SignatureSets) 
- **Access Method**: Direct download or programmatic retrieval via API (if applicable)
- **Data Format**: rda, CSV (signatures, metadata)

Signatures were applied using rank-based, enrichment-based, or weighted scoring approaches, depending on the analysis.

---

## Key Clinical Variables

Clinical annotations were harmonized across datasets when applicable. Key variables include:

| Variable            | Description                                       | Format   | Example      |
|---------------------|---------------------------------------------------|----------|--------------|
| patientid           | Unique patient identifier                         | string   | GSE12345_P01 |
| age                 | Age at diagnosis                                  | integer  | 63           |
| sex                 | Biological sex                                    | factor   | M/F          |
| cancer_type         | Primary cancer type                               | string   | Melanoma     |
| histo               | Histological classification                       | string   | Melanoma     |
| treatment_type      | IO therapy category                               | string   | PD-1/PD-L1   |
| stage               | Tumor stage at diagnosis                          | string   | Stage II     |
| recist              | RECIST clinical response                          | factor   | CR/PR/SD/PD  |
| response            | Clinical benefit status (e.g., response)          | string   | R/NR         |
| survival_time_os    | Overall survival time (months)                    | numeric  | 21.3         |
| event_occurred_os   | Overall survival event (1 = death)                | binary   | 1            |
| survival_time_pfs   | Progression-free survival time (months)           | numeric  | 18.2         |
| event_occurred_pfs  | Progression-free survival event (1 = progression) | binary   | 1            |
| survival_unit       | Unit of survival time                             | string   | months       |

---
