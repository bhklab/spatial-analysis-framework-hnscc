suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(clustree)
})

# ============================================================
# Part 2.1 — Non-malignant spot characterization
# - Refine non-malignant compartments by reclustering:
#   Fibroblasts, Macrophages, B cells 
# - Keeps malignant labels as-is (Cancer_C*)
# - Does NOT assign biological names; keeps *_C0, *_C1, ...
# ============================================================

# ---- EDIT PATHS ----
in_rds  <- "/part2/01_malignant_characterization/RM_ST_with_malignant_subclusters.rds"
out_dir <- "/part2/02_nonmalignant_characterization"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- PARAMETERS ----
dims_use <- 1:15

# choose resolutions (tweak later; these are safe starting points)
fibro_res <- 0.12
macro_res <- 0.12
bcell_res <- 0.20

# resolution scan for QC with clustree
res_scan <- c(seq(0.02, 0.40, by = 0.02), seq(0.50, 1.0, by = 0.10))

# -----------------------------
# Helpers
# -----------------------------
scan_resolutions <- function(obj, res_vec) {
  tmp <- obj
  for (r in res_vec) tmp <- FindClusters(tmp, resolution = r, verbose = FALSE)
  tmp
}

label_as_C <- function(prefix, x) {
  paste0(prefix, "_C", as.integer(as.character(x)))
}

# -----------------------------
# Load
# -----------------------------
RM_ST <- readRDS(in_rds)

# If layers exist, join (safe)
RM_ST <- tryCatch(JoinLayers(RM_ST), error = function(e) RM_ST)

# Use SCT for clustering/UMAP if available
if (!"SCT" %in% names(RM_ST@assays)) {
  stop("RM_ST is missing SCT assay. Run Part 2.0 (01a/01b) first.")
}
DefaultAssay(RM_ST) <- "SCT"

# This column should exist from Part 2.0
stopifnot("Final_Malignant_SubCluster" %in% colnames(RM_ST@meta.data))

# Initialize refined label column
RM_ST$Spatial_CT_refined <- as.character(RM_ST$Final_Malignant_SubCluster)

# -----------------------------
# 2.1A — Fibroblast refinement
# -----------------------------
Fibro <- subset(RM_ST, subset = Spatial_CT_refined %in% "Fibroblasts")
if (ncol(Fibro) > 0) {
  
  # require integrated embedding from Part 2.0
  if (!"integrated.dr" %in% names(Fibro@reductions)) {
    stop("Fibro subset missing integrated.dr reduction. Ensure 01a ran successfully.")
  }
  
  Fibro <- RunUMAP(
    Fibro, reduction = "integrated.dr", dims = dims_use,
    reduction.name = "umap.fibro", verbose = FALSE
  )
  Fibro <- FindNeighbors(Fibro, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)
  Fibro <- FindClusters(Fibro, resolution = fibro_res, cluster.name = "Fibro_clusters", verbose = FALSE)
  
  # scan for QC
  Fibro_scan <- scan_resolutions(Fibro, res_scan)
  saveRDS(Fibro_scan, file.path(out_dir, "Fibro_resolution_scan.rds"))
  
  # label: Fibro_C0, Fibro_C1, ...
  Fibro$Fibro_refined <- label_as_C("Fibro", Fibro$Fibro_clusters)
  
  # push back to RM_ST
  common <- intersect(colnames(RM_ST), colnames(Fibro))
  RM_ST$Spatial_CT_refined[common] <- Fibro$Fibro_refined[common]
  
  saveRDS(Fibro, file.path(out_dir, "Fibro_refined.rds"))
}

# -----------------------------
# 2.1B — Macrophage refinement
# -----------------------------
Macro <- subset(RM_ST, subset = Spatial_CT_refined %in% "Macrophages")
if (ncol(Macro) > 0) {
  
  if (!"integrated.dr" %in% names(Macro@reductions)) {
    stop("Macrophage subset missing integrated.dr reduction. Ensure 01a ran successfully.")
  }
  
  Macro <- RunUMAP(
    Macro, reduction = "integrated.dr", dims = dims_use,
    reduction.name = "umap.macro", verbose = FALSE
  )
  Macro <- FindNeighbors(Macro, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)
  Macro <- FindClusters(Macro, resolution = macro_res, cluster.name = "Macro_clusters", verbose = FALSE)
  
  Macro_scan <- scan_resolutions(Macro, res_scan)
  saveRDS(Macro_scan, file.path(out_dir, "Macro_resolution_scan.rds"))
  
  Macro$Macro_refined <- label_as_C("Macro", Macro$Macro_clusters)
  
  common <- intersect(colnames(RM_ST), colnames(Macro))
  RM_ST$Spatial_CT_refined[common] <- Macro$Macro_refined[common]
  
  saveRDS(Macro, file.path(out_dir, "Macro_refined.rds"))
}

# -----------------------------
# 2.1C — B cell refinement 
# -----------------------------
B <- subset(RM_ST, subset = Spatial_CT_refined %in% "B cells")
if (ncol(B) > 0) {
  
  if (!"integrated.dr" %in% names(B@reductions)) {
    stop("B cell subset missing integrated.dr reduction. Ensure 01a ran successfully.")
  }
  
  B <- RunUMAP(
    B, reduction = "integrated.dr", dims = dims_use,
    reduction.name = "umap.bcells", verbose = FALSE
  )
  B <- FindNeighbors(B, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)
  B <- FindClusters(B, resolution = bcell_res, cluster.name = "B_clusters", verbose = FALSE)
  
  B_scan <- scan_resolutions(B, res_scan)
  saveRDS(B_scan, file.path(out_dir, "Bcell_resolution_scan.rds"))
  
  B$B_refined <- label_as_C("B", B$B_clusters)
  
  common <- intersect(colnames(RM_ST), colnames(B))
  RM_ST$Spatial_CT_refined[common] <- B$B_refined[common]
  
  saveRDS(B, file.path(out_dir, "Bcell_refined.rds"))
}

# -----------------------------
# Save updated RM_ST
# -----------------------------
saveRDS(RM_ST, file.path(out_dir, "RM_ST_with_nonmalignant_refined_labels.rds"))
message("Saved: ", file.path(out_dir, "RM_ST_with_nonmalignant_refined_labels.rds"))
