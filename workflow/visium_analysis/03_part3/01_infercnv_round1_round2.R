#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(infercnv)
})

# ============================================================
# Part 3 / inferCNV — Round 1 + Round 2
#   Round 1: malignant vs non_malignant
#   Round 2: malignant subclusters vs non_malignant
#
# Input:
#   - Seurat object (RM_ST) containing:
#       * Broad_CT (or malignancy flag)
#       * malignant subcluster column for Round 2 (e.g., Final_Malignant_SubCluster)
# Output:
#   - outputs/part3/infercnv/round1/...
#   - outputs/part3/infercnv/round2/...
# ============================================================

# -------------------------
# EDIT PATHS
# -------------------------
in_rds   <- "part3/RM_ST.rds"                       # from Part 2
out_root <- "part3/infercnv"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# inferCNV requires a gene order file: gene \t chr \t start \t end
gene_order_file <- "path/to/gene_order.tsv"
if (!file.exists(gene_order_file)) {
  stop("Missing gene_order.tsv: ", gene_order_file,
       "\nCreate a 4-col file: gene, chr, start, end.")
}

# -------------------------
# PARAMETERS
# -------------------------
assay_use   <- "Spatial"
cutoff      <- 0.1
num_threads <- 8

# How to define malignancy for Round 1
# (This is the most common rule in your pipeline)
broad_ct_col <- "Broad_CT"
malignant_value <- "Malignant cells"

# Malignant subcluster column for Round 2
# (Keep IDs like Cancer_C0, Cancer_C1... for now; no biological naming)
malig_subcluster_col <- "Final_Malignant_SubCluster"

# Reference group name (must exist in annotations)
ref_group <- "non_malignant"

# -------------------------
# Helpers
# -------------------------
write_annotations <- function(obj, groups, out_file) {
  stopifnot(length(groups) == ncol(obj))
  anno <- data.frame(cell_id = colnames(obj), group = groups, stringsAsFactors = FALSE)
  write.table(anno, out_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  invisible(out_file)
}

run_infercnv_one <- function(counts_mat, anno_file, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_mat,
    annotations_file  = anno_file,
    delim             = "\t",
    gene_order_file   = gene_order_file,
    ref_group_names   = ref_group
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff            = cutoff,
    out_dir           = out_dir,
    cluster_by_groups = TRUE,
    denoise           = TRUE,
    HMM               = TRUE,
    num_threads       = num_threads
  )
  
  infercnv_obj
}

# -------------------------
# Load object
# -------------------------
RM_ST <- readRDS(in_rds)
DefaultAssay(RM_ST) <- assay_use

md <- RM_ST@meta.data
stopifnot(broad_ct_col %in% colnames(md))

counts_mat <- GetAssayData(RM_ST, assay = assay_use, slot = "counts")

# ============================================================
# ROUND 1: malignant vs non_malignant
# ============================================================
round1_dir <- file.path(out_root, "round1")
dir.create(round1_dir, recursive = TRUE, showWarnings = FALSE)

groups_round1 <- ifelse(md[[broad_ct_col]] == malignant_value, "malignant", "non_malignant")
anno1_file <- file.path(round1_dir, "annotations_round1.txt")
write_annotations(RM_ST, groups_round1, anno1_file)

out1_run <- file.path(round1_dir, "infercnv_run_round1")
infercnv_obj_r1 <- run_infercnv_one(counts_mat, anno1_file, out1_run)
saveRDS(infercnv_obj_r1, file.path(round1_dir, "infercnv_round1_obj.rds"))

message("✅ Round 1 done: ", out1_run)

# ============================================================
# ROUND 2: malignant subclusters vs non_malignant
#   - malignant cells split by malig_subcluster_col
#   - non-malignant stays "non_malignant" (reference)
# ============================================================
round2_dir <- file.path(out_root, "round2")
dir.create(round2_dir, recursive = TRUE, showWarnings = FALSE)

if (!(malig_subcluster_col %in% colnames(md))) {
  stop("Missing malignant subcluster column for Round 2: ", malig_subcluster_col,
       "\nAdd it in Part 2 first (keep Cancer_C* IDs), then rerun.")
}

groups_round2 <- md[[broad_ct_col]]
groups_round2 <- ifelse(
  md[[broad_ct_col]] == malignant_value,
  as.character(md[[malig_subcluster_col]]),  # e.g., Cancer_C0, Cancer_C1...
  ref_group
)

# Safety: no NAs for inferCNV
if (any(is.na(groups_round2)) || any(groups_round2 == "")) {
  stop("Round 2 annotations contain NA/empty values. Check ", malig_subcluster_col)
}

anno2_file <- file.path(round2_dir, "annotations_round2.txt")
write_annotations(RM_ST, groups_round2, anno2_file)

out2_run <- file.path(round2_dir, "infercnv_run_round2")
infercnv_obj_r2 <- run_infercnv_one(counts_mat, anno2_file, out2_run)
saveRDS(infercnv_obj_r2, file.path(round2_dir, "infercnv_round2_obj.rds"))

message("✅ Round 2 done: ", out2_run)

message("🎉 Finished inferCNV Round 1 + Round 2. Outputs in: ", out_root)
