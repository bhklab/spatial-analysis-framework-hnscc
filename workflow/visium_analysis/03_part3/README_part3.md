
# Part 3 — CNV inference (inferCNV) and cell–cell communication (CellChat)

This folder contains the **Part 3** analysis code used for the manuscript:
1) **inferCNV** (two-round CNV inference)
2) **CellChat** (multi-tier communication inference; compute-only)

**Important:** This part is intentionally **compute + saved outputs only**.  
All **visualization/figure generation** is handled in **Part 4**.

---

## Expected inputs (from Part 1–2)

You should already have:

- `RM_ST.rds`  
  A Seurat object containing all spots with finalized metadata from Part 2 (e.g., `Broad_CT`, `Spatial_CT`, etc.)

- `Cancer_v3.RDS` (or equivalent malignant-only Seurat object)  
  A Seurat object restricted to malignant spots for malignant-only CellChat.

- `gene_order.tsv` (required by inferCNV)  
  A 4-column TSV with:
  - `gene` (symbol matching the expression matrix rownames)
  - `chr`
  - `start`
  - `end`

---

## Outputs

### inferCNV
Written under:
- `outputs/part3/infercnv/round1/`
- `outputs/part3/infercnv/round2/`

Each round contains:
- `annotations_round*.txt` (cell/spot → group labels for inferCNV)
- `infercnv_run_round*/` (inferCNV run directory with standard outputs)
- `infercnv_round*_obj.rds` (saved inferCNV object for downstream use)

Round definitions:
- **Round 1:** `malignant` vs `non_malignant`
- **Round 2:** malignant split into **subclusters** vs `non_malignant`

---

### CellChat
Written under:
- `outputs/part3/cellchat/tier1_malignant_only/`
- `outputs/part3/cellchat/tier2_global/`
- `outputs/part3/cellchat/tier2b_malignant_plus_caf/`
- `outputs/part3/cellchat/tier2c_malignant_plus_immune/`

Each tier contains:
- `cellchat.rds` (full CellChat object)
- `communications_all.csv` (ligand–receptor interactions)
- `communications_pathway_level.csv` (pathway-level interactions; if available)
- `net_count_matrix.csv` (interaction counts)
- `net_weight_matrix.csv` (interaction weights)

**Note:** No plots are generated here by design (Part 4 handles plotting).

---

## Script inventory and execution order

### 1) inferCNV (Round 1 + Round 2)
Run:
- `code/03_part3/01_infercnv/01a_infercnv_round1_round2.R`

This script:
- loads `RM_ST.rds`
- creates inferCNV annotations
- runs inferCNV twice:
  - Round 1: malignant vs non-malignant
  - Round 2: malignant subclusters vs non-malignant
- saves inferCNV objects and run directories for both rounds

**You must set these paths inside the script:**
- `in_rds`
- `out_root`
- `gene_order_file`

---

### 2) CellChat (Tier1 + Tier2 + Tier2b + Tier2c)
Run:
- `code/03_part3/02_cellchat/01_run_cellchat_tier1_tier2_tier2b_tier2c.R`

This script:
- loads `Cancer_v3.RDS` (Tier1) and `RM_ST.rds` (Tier2/Tier2b/Tier2c)
- runs CellChat pipelines in four modes:
  - **Tier1:** malignant-only (grouped by malignant cluster IDs)
  - **Tier2:** global TME model (uses a relabeled/condensed grouping column)
  - **Tier2b:** malignant + CAF only
  - **Tier2c:** malignant + immune only
- saves CellChat objects and key result tables for downstream visualization

**You must set these paths inside the script:**
- `rm_st_rds`
- `cancer_rds`
- `out_root`

---

## Notes on labels / cluster naming

- For inferCNV Round 2 and CellChat Tier1/Tier2, we intentionally keep **cluster IDs unchanged** (e.g., `Cancer_C0`, `Cancer_C1`, …).
- Biological naming/interpretation of clusters is deferred to **Part 4** to keep computation reproducible and avoid figure-driven renaming in core analysis code.

---

## Reproducibility tips

- Avoid `setwd()` in scripts. Use explicit input/output paths.
- Save intermediate objects (`*.rds`) after major compute steps.
- Keep `gene_order.tsv` versioned or documented (genome build, source).
- If running on HPC, use job scripts to run inferCNV/CellChat and store logs alongside outputs.

---

## Troubleshooting

- **inferCNV fails with missing genes:** ensure `gene_order.tsv` gene symbols match the expression matrix rownames.
- **CellChat drops groups:** increase/decrease `min_cells_per_group` in the CellChat script.
- **Large memory usage:** consider DietSeurat before CellChat, or run tiers separately on HPC.
