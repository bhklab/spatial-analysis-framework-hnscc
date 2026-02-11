suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

# ============================================================
# Part 2.0 (01a) — Global SCT + RPCA integration + UMAPs
# Input:  Part 1 output (AllST.rds or RM_ST.rds equivalent)
# Output: RM_ST_integrated.rds (with umap.unintegrated + umap.integrated)
# ============================================================

# ---- EDIT PATHS ----
in_rds  <- "/part1/08_final_malignant_annotation/AllST.rds"
out_dir <- "/part2/01_malignant_characterization"
out_rds <- file.path(out_dir, "RM_ST_integrated.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- PARAMETERS ----
dims_use <- 1:15
k_weight <- 50
plan("sequential")  # safe default; change if you use multicore
options(future.globals.maxSize = 30 * 1024^3)

# -----------------------------
# Load
# -----------------------------
RM_ST <- readRDS(in_rds)

# -----------------------------
# Ensure Spatial assay is ready
# -----------------------------
DefaultAssay(RM_ST) <- "Spatial"

# -----------------------------
# SCT + PCA + unintegrated UMAP
# -----------------------------
RM_ST <- SCTransform(RM_ST, assay = "Spatial", verbose = FALSE)
RM_ST <- RunPCA(RM_ST, verbose = FALSE)

RM_ST <- RunUMAP(
  RM_ST,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap.unintegrated",
  verbose = FALSE
)

RM_ST <- FindNeighbors(RM_ST, reduction = "pca", dims = dims_use, verbose = FALSE)

# Light clustering just for QC visualization 
RM_ST <- FindClusters(
  RM_ST,
  resolution = 0.05,
  cluster.name = "unintegrated_clusters",
  verbose = FALSE
)

# -----------------------------
# RPCA integration (Seurat v5 IntegrateLayers)
# -----------------------------
RM_ST <- IntegrateLayers(
  object = RM_ST,
  method = RPCAIntegration,
  normalization.method = "SCT",
  k.weight = k_weight
)

# Neighbor graph + integrated UMAP
RM_ST <- FindNeighbors(RM_ST, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)

RM_ST <- RunUMAP(
  RM_ST,
  reduction = "integrated.dr",
  dims = dims_use,
  reduction.name = "umap.integrated",
  verbose = FALSE
)

# -----------------------------
# Save
# -----------------------------
saveRDS(RM_ST, out_rds)
message("Saved: ", out_rds)
