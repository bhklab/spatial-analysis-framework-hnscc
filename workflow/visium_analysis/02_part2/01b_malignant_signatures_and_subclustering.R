suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(AUCell)
  library(readxl)
  library(clustree)
})

# ============================================================
# Part 2.0 (01b) — Malignant signatures + multi-stage subclustering
# Input:  RM_ST_integrated.rds from 01a
# Output: Cancer objects + RM_ST_with_malignant_subclusters.rds
# Notes:  Does NOT assign biological names; keeps Cancer_C0 etc.
# ============================================================

# ---- EDIT PATHS ----
in_rds   <- "/part2/01_malignant_characterization/RM_ST_integrated.rds"
deg_xlsx <- "/DEGs_BoseLabSignature.xlsx"
out_dir  <- "/part2/01_malignant_characterization"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- PARAMETERS ----
dims_use <- 1:15

# multi-stage parameters (match your script intent; tweak later if needed)
lvl1_res <- 0.05
lvl2_res <- 0.05
lvl3_res <- 0.23

# For clustree QA scans (saved as rds)
res_scan <- c(seq(0.01, 0.40, by = 0.01), seq(0.50, 1.0, by = 0.10))

# Which malignant subset column to use
malignant_label_col <- "Broad_CT"
malignant_label_val <- "Malignant cells"

# -----------------------------
# Load
# -----------------------------
RM_ST <- readRDS(in_rds)

stopifnot(malignant_label_col %in% colnames(RM_ST@meta.data))

# If layers exist, join (safe)
RM_ST <- tryCatch(JoinLayers(RM_ST), error = function(e) RM_ST)

# -----------------------------
# AUCell signature scoring (adds AUC_* to metadata)
# -----------------------------
DEGs <- readxl::read_excel(deg_xlsx)

core_genes <- DEGs %>% filter(group == "n1") %>% pull(gene)
edge_genes <- DEGs %>% filter(group == "n3") %>% pull(gene)

core_genes <- intersect(core_genes, rownames(RM_ST))
edge_genes <- intersect(edge_genes, rownames(RM_ST))

Basal_markers <- c("KRT5","KRT14","TP63","ITGA6","DSG3","LAMB3","LAMC2","S100A2")
Classical_markers <- c("KRT8","KRT18","KRT19","SFN","CLDN4","EPCAM","LY6D","AGR2","MUC1","GRHL3","TJP3")
Differentiated_markers <- c("IVL","SPRR1A","SPRR1B","SPRR2A","SPRR2D","TGM1","TGM3","KRT1","KRT10","LOR","FLG")
EMT_markers <- c("VIM","FN1","SNAI2","TWIST1","ZEB1","CDH2","COL1A1","COL1A2","COL3A1","COL5A1")
Proliferative_markers <- c("MKI67","PCNA","TOP2A","MCM2","MCM4","AURKA","BIRC5")
Hypoxic_markers <- c("CA9","VEGFA","LDHA","HIF1A","ENO2","NDRG1")
Inflammatory_markers <- c("CXCL9","CXCL10","CXCL11","STAT1","IFI6")

gene_sets <- list(
  Tumor_Core    = core_genes,
  Leading_Edge  = edge_genes,
  Basal         = Basal_markers,
  Classical     = Classical_markers,
  Differentiated= Differentiated_markers,
  EMT           = EMT_markers,
  Proliferative = Proliferative_markers,
  Hypoxic       = Hypoxic_markers,
  Inflammatory  = Inflammatory_markers
)

# AUCell expects genes x cells matrix.
DefaultAssay(RM_ST) <- "Spatial"
counts <- GetAssayData(RM_ST, assay = "Spatial", slot = "counts")

cell_sums <- Matrix::colSums(counts)
counts_cpm <- t(t(counts) / cell_sums * 10000)
exprMatrix <- as(counts_cpm, "dgCMatrix")  # genes x spots

rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE, verbose = FALSE)
auc_obj  <- AUCell_calcAUC(gene_sets, rankings, verbose = FALSE)

auc_mat <- t(getAUC(auc_obj))  # spots x geneSets
auc_df <- as.data.frame(auc_mat)
colnames(auc_df) <- paste0("AUC_", colnames(auc_df))

RM_ST <- AddMetaData(RM_ST, metadata = auc_df)

# -----------------------------
# Cancer subset
# -----------------------------
Cancer <- subset(RM_ST, subset = .data[[malignant_label_col]] %in% malignant_label_val)
if (ncol(Cancer) == 0) stop("No malignant spots found using: ", malignant_label_col, " == ", malignant_label_val)

DefaultAssay(Cancer) <- "SCT"

# ensure we have the integrated reduction
if (!"integrated.dr" %in% names(Cancer@reductions)) {
  stop("Cancer object lacks integrated.dr reduction. Run 01a first or ensure IntegrateLayers worked.")
}

# Stage-0 UMAP for Cancer (for visual checks)
Cancer <- RunUMAP(
  Cancer,
  reduction = "integrated.dr",
  dims = dims_use,
  reduction.name = "umap.cancer_lvl1",
  verbose = FALSE
)
Cancer <- FindNeighbors(Cancer, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)

# -----------------------------
# Helper: resolution scan + clustree object
# -----------------------------
scan_resolutions <- function(obj, res_vec, prefix = "SCT_snn_res.") {
  tmp <- obj
  for (r in res_vec) tmp <- FindClusters(tmp, resolution = r, verbose = FALSE)
  # clustree uses prefix matching default column naming from FindClusters
  tmp
}

# -----------------------------
# Level 1 clustering
# -----------------------------
Cancer <- FindClusters(Cancer, resolution = lvl1_res, cluster.name = "Cancer_lvl1", verbose = FALSE)

# scan object for lvl1 QA
Cancer_lvl1_scan <- scan_resolutions(Cancer, res_scan)
saveRDS(Cancer_lvl1_scan, file.path(out_dir, "Cancer_lvl1_resolution_scan.rds"))

# -----------------------------
# Level 2: subset lvl1 cluster "0" and recluster
# -----------------------------
Cancer_lvl2 <- subset(Cancer, subset = Cancer_lvl1 %in% "0")

Cancer_lvl2 <- RunUMAP(
  Cancer_lvl2,
  reduction = "integrated.dr",
  dims = dims_use,
  reduction.name = "umap.cancer_lvl2",
  verbose = FALSE
)
Cancer_lvl2 <- FindNeighbors(Cancer_lvl2, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)
Cancer_lvl2 <- FindClusters(Cancer_lvl2, resolution = lvl2_res, cluster.name = "Cancer_lvl2", verbose = FALSE)

Cancer_lvl2_scan <- scan_resolutions(Cancer_lvl2, res_scan)
saveRDS(Cancer_lvl2_scan, file.path(out_dir, "Cancer_lvl2_resolution_scan.rds"))

# -----------------------------
# Level 3: subset lvl2 cluster "0" and recluster at higher res
# -----------------------------
Cancer_lvl3 <- subset(Cancer_lvl2, subset = Cancer_lvl2 %in% "0")

Cancer_lvl3 <- RunUMAP(
  Cancer_lvl3,
  reduction = "integrated.dr",
  dims = dims_use,
  reduction.name = "umap.cancer_lvl3",
  verbose = FALSE
)
Cancer_lvl3 <- FindNeighbors(Cancer_lvl3, reduction = "integrated.dr", dims = dims_use, verbose = FALSE)
Cancer_lvl3 <- FindClusters(Cancer_lvl3, resolution = lvl3_res, cluster.name = "Cancer_lvl3", verbose = FALSE)

Cancer_lvl3_scan <- scan_resolutions(Cancer_lvl3, res_scan)
saveRDS(Cancer_lvl3_scan, file.path(out_dir, "Cancer_lvl3_resolution_scan.rds"))

# -----------------------------
# Propagate labels upward WITHOUT biological names
# Final label column: Final_Malignant_SubCluster
# - Start from lvl3 clusters ("0","1",...) -> make "Cancer_C0","Cancer_C1",...
# - Push into lvl2 -> push into lvl1 -> push into Cancer
# -----------------------------
make_cancer_c_labels <- function(x) paste0("Cancer_C", as.integer(as.character(x)))

Cancer_lvl3$Final_Malignant_SubCluster <- make_cancer_c_labels(Cancer_lvl3$Cancer_lvl3)

# lvl2 default label = Cancer_C<lvl2>
Cancer_lvl2$Final_Malignant_SubCluster <- make_cancer_c_labels(Cancer_lvl2$Cancer_lvl2)
common_23 <- intersect(colnames(Cancer_lvl2), colnames(Cancer_lvl3))
Cancer_lvl2$Final_Malignant_SubCluster[common_23] <- Cancer_lvl3$Final_Malignant_SubCluster[common_23]

# lvl1 default label = Cancer_C<lvl1>
Cancer$Final_Malignant_SubCluster <- make_cancer_c_labels(Cancer$Cancer_lvl1)
common_12 <- intersect(colnames(Cancer), colnames(Cancer_lvl2))
Cancer$Final_Malignant_SubCluster[common_12] <- Cancer_lvl2$Final_Malignant_SubCluster[common_12]

# -----------------------------
# Push malignant subclusters back into RM_ST
# (non-malignant stays as Broad_CT or "Other"; we keep your original behavior)
# -----------------------------
RM_ST$Final_Malignant_SubCluster <- as.character(RM_ST[[malignant_label_col]][,1])
common_rm <- intersect(colnames(RM_ST), colnames(Cancer))
RM_ST$Final_Malignant_SubCluster[common_rm] <- Cancer$Final_Malignant_SubCluster[common_rm]

# -----------------------------
# Save outputs
# -----------------------------
saveRDS(Cancer,      file.path(out_dir, "Cancer_lvl1_with_final_subclusters.rds"))
saveRDS(Cancer_lvl2, file.path(out_dir, "Cancer_lvl2_refined.rds"))
saveRDS(Cancer_lvl3, file.path(out_dir, "Cancer_lvl3_refined.rds"))
saveRDS(RM_ST,       file.path(out_dir, "RM_ST_with_malignant_subclusters.rds"))

message("Done.")
message("Saved Cancer lvl1: ", file.path(out_dir, "Cancer_lvl1_with_final_subclusters.rds"))
message("Saved RM_ST:      ", file.path(out_dir, "RM_ST_with_malignant_subclusters.rds"))
