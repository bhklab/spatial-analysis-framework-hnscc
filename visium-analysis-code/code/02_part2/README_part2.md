
# Part 2 — Spot Characterization, Neighborhoods, and Near vs Far DE

This section performs:
1) global integration and embedding of the spatial dataset,  
2) malignant spot refinement via multi-stage subclustering (**no biological naming yet**),  
3) non-malignant compartment refinement (**no biological naming yet**),  
4) neighborhood analysis around malignant states,  
5) defining *near* vs *far* spots for each (neighbor cell type N, malignant state M) pair,  
6) differential expression (DE) in N: **near vs far**, conditioned on malignant state M.

**Note:** Biological interpretation / naming of clusters is deferred to **Part 4 (Visualization & annotation)**.  
Clusters are kept as `Cancer_C*`, `Fibro_C*`, `Macro_C*`, etc.

---

## Inputs

### Required
- **Part 1 output Seurat object** (RDS): a combined spatial object with metadata (e.g., `ID`, `Broad_CT`) and Visium spatial coordinates.

### Additional (for AUCell scoring in 2.0)
- `DEGs_BoseLabSignature.xlsx`
  - Must contain at least two groups:
    - `n1` (Tumor Core)
    - `n3` (Leading Edge)

---

## Outputs

### 2.0 Malignant characterization (integration + subclustering)
- `RM_ST_integrated.rds`
- `RM_ST_with_malignant_subclusters.rds`
- `Cancer_lvl1_with_final_subclusters.rds`
- `Cancer_lvl2_refined.rds`
- `Cancer_lvl3_refined.rds`

### 2.1 Non-malignant characterization
- `RM_ST_with_nonmalignant_refined_labels.rds`

### 2.2 Neighborhood analysis
- `combined_neighborhood_counts.csv`
- `Sample_Summary_Table.csv`

### 2.3 Identify near vs far (graph-distance based)
- `annotated_objects_list_forDEGs_withFar.rds`
  - Adds `*_far_*` columns corresponding to each `*_neighbouring_*` pair.

### 2.4 Differential expression: near vs far
- `DEG_near_vs_far/DEG_<sample>.csv` (per sample, per pair)
- `DEG_near_vs_far/OVERLAP_<pair>.csv`
- `DEG_near_vs_far/OVERLAP_all_pairs_summary.csv`
- `DEG_near_vs_far/DEG_near_vs_far_results.rds` (bundle of results + parameters)

---

## Script order

Run scripts in this order:

### 2.0 — Malignant spot characterization
1. `code/02_part2/01a_global_integration_and_umap.R`
2. `code/02_part2/01b_malignant_signatures_and_subclustering.R`

### 2.1 — Non-malignant spot characterization
3. `code/02_part2/02_nonmalignant_spot_characterization.R`

### 2.2 — Neighborhood analysis
4. `code/02_part2/03_neighborhood_analysis.R`

### 2.3 — Identify far spots for each (N, M) pair
5. `code/02_part2/04_identifying_near_far.R`

### 2.4 — Differential expression (near vs far)
6. `code/02_part2/05_dge_near_vs_far.R`

> **Note:** Negative control scripts (e.g., prior `2.2b`) are intentionally excluded from this publication pipeline.

---

## Key analysis notes

### Multi-stage malignant refinement (Part 2.0)
Malignant refinement is performed in stages:
- **Level 1:** cluster all malignant spots (`Cancer_lvl1`)
- **Level 2:** subset `Cancer_lvl1 == "0"` and recluster (`Cancer_lvl2`)
- **Level 3:** subset `Cancer_lvl2 == "0"` and recluster (`Cancer_lvl3`)

Final labels are propagated upward and stored as:
- `Final_Malignant_SubCluster` with values `Cancer_C0`, `Cancer_C1`, ...

### Non-malignant refinement (Part 2.1)
Selected non-malignant compartments are refined independently:
- Fibroblasts → `Fibro_C*`
- Macrophages → `Macro_C*`
- B cells (optional) → `B_C*`

Result is stored in:
- `Spatial_CT_refined` (or equivalent label column used downstream)

### Neighborhood analysis (Part 2.2)
Neighborhoods are computed using **k-nearest neighbors (k=6)** using spatial spot coordinates.
Neighborhood counting focuses on:
- malignant states (`Cancer_C*`) as the “primary”
- non-malignant labels as “neighbors”

This step typically produces `*_neighbouring_*` metadata columns used downstream.

### Far definition (Part 2.3)
For each pair (neighbor cell type **N**, malignant state **M**):
- **Near** spots are N spots labeled as neighboring M (`N_neighbouring_M == 1`)
- **Far** spots are N spots that are:
  - not near M, and
  - sufficiently distant from M in the neighbor graph (minimum hop threshold), with optional quantile fallback.

This creates `N_far_M` columns aligned with the existing `N_neighbouring_M` columns.

### Differential expression: near vs far (Part 2.4)
Within each sample, for each (N, M) pair:
- Compare N spots that are **near M** vs **far from M**.
- DE is performed per sample and exported as CSV.
- Overlap summaries across samples are exported per pair and as a combined table.

---

## Expected metadata fields

These scripts assume the following are present in the Seurat object metadata:
- `ID` (sample identifier used for `SplitObject`)
- `Broad_CT` (broad cell-type assignment including `Malignant cells`)
- One label column used consistently downstream (e.g., `Spatial_CT` or `Spatial_CT_refined`)
- Spatial coordinates available via the object’s `images` slot (Visium)

If your object uses different column names, update the corresponding parameters at the top of each script.

---

## Reproducibility

- Do **not** commit large processed RDS objects to GitHub.
- Store processed objects and large outputs in local/project storage.
- Keep scripts deterministic and parameterized via editable paths at the top of each script.
- Export tables (CSV/TSV) that allow reviewers to regenerate figures.
