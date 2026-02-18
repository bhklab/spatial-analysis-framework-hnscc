# Scripts Directory - External IO Validation RNAseq

## Purpose

This directory provides **modular R scripts** implementing a **spatially informed gene signature analysis pipeline** for immunotherapy transcriptomic datasets.

The workflow integrates:

- Rank-based single-sample gene signature scoring  
- Centroid-based spatial projection  
- Association testing with clinical outcomes  
- Cross-cohort meta-analysis  
- Signature correlation and clustering  
- Publication-ready visualization  

All analyses are modular and reproducible.

---

## Key Capabilities

- Single-sample signature scoring (singscore, GSVA, ssGSEA, weighted mean)
- Centroid-distance–based spatial signature scoring
- Cox proportional hazards and logistic regression analyses
- Pan-cancer and cancer-specific meta-analysis
- Signature correlation and meta-correlation
- Gene-overlap–based clustering with pathway enrichment
- Automated generation of publication-ready figures

---

# Modular Scripts (`workflow/external_validation_rnaseq/`)

The workflow is organized into eight modular components:

---

## (1) `runProcData.R`

Prepares transcriptomic datasets for downstream analysis.

**Tasks:**
- Load `MultiAssayExperiment` objects  
- Convert to `SummarizedExperiment`  
- Normalize expression matrices  
- Harmonize clinical variables  
- Load and validate curated signature libraries  

---

## (2) `runSigScore.R`

Computes gene signature scores.

**Implemented Methods:**
- `singscore`
- GSVA
- ssGSEA
- Weighted mean
- Specific methods (e.g., IPS, COX-IS)
- Spatial centroid-based scoring

### Spatially Informed Scoring

For spatial signatures:

1. Perform rank-based single-sample scoring  
2. Project samples into a centroid-defined signature space  
3. Compute Euclidean distances to two centroids  
4. Define the final score as the distance to **Centroid 2**

Centroid 2 represents a spatial tumor microenvironment state associated with favorable progression-free survival under immune checkpoint blockade (ICB).

Lower distance indicates greater similarity to the favorable spatial state.

---

## (3) `runSigAssoc.R`

Performs association analyses between signature scores and clinical outcomes.

**Analyses:**
- Cox proportional hazards models (PFS)
- Logistic regression models (response)
- Continuous predictors
- Median-based dichotomized predictors

---

## (4) `runMeta.R`

Performs cross-study meta-analysis of association results.

**Features:**
- Fixed-effects and random-effects models  
- Pan-cancer integration  
- Cancer-specific analysis  
- Treatment-specific analysis  
- Heterogeneity metrics (I², Q-test)  

**Requirement:** ≥3 independent cohorts per feature.

---

## (5) `runCorr.R`

Performs correlation and meta-correlation analysis of signature scores.

**Includes:**
- Pearson correlations per dataset  
- Meta-correlation across studies  
- Hierarchical clustering  
- Heatmap generation  

**Requirement:** ≥3 datasets per signature pair.

---

## (6) `runSigCluster.R`

Clusters gene signatures based on shared-gene overlap.

**Workflow:**
1. Build signature–gene overlap matrix  
2. Perform PCA  
3. Apply affinity propagation clustering  
4. Conduct KEGG enrichment analysis (Enrichr)  

Includes overlap matrices, PCA plots, clustering heatmaps, and pathway enrichment tables.

---

## (7) `sigDistanceFunction.R`

Implements centroid-based spatial distance scoring.

**Functions:**
- `compute_distances()` — Computes Euclidean distances to predefined centroids  
- `geneSigNovel()` — Performs rank-based scoring and returns distance to Centroid 2  

**Supported Input Formats:**
- `SummarizedExperiment`
- `MultiAssayExperiment`
- matrix/data.frame (genes × samples)

---

## (8) `runVisualization.R`

Generates publication-ready figures.

**Figures include:**
- Kaplan–Meier survival curves  
- Boxplots (Responder vs Non-responder)  
- Forest plots (logHR and logOR)  
- Meta-analysis summary visualizations  


---

# Input Requirements

Each dataset should provide:

- A normalized gene expression matrix (e.g., log2(TPM + 1))  
- Clinical variables:
  - Progression-Free Survival (time and event)
  - Binary response status  
- Curated gene signature definitions  
- Predefined centroid profiles (.rds files) for spatial signatures  

---

# Notes

- PFS is administratively censored at 24 months.  
- Median dichotomization requires ≥3 samples per group.  
- Meta-analysis requires ≥3 independent cohorts.  
- Correlation meta-analysis requires ≥3 datasets per signature pair.  

---

