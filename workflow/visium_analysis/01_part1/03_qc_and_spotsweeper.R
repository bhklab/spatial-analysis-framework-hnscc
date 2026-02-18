suppressPackageStartupMessages({
  library(Seurat)
  library(SpotSweeper)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(dplyr)
})

# ============================================================
# Part 1.2 — QC + SpotSweeper local outliers + global cutoffs
# Input:
#   - in_dir:  folder containing *.rds Seurat objects (from Part 1.1)
# Output:
#   - out_dir: folder containing filtered *.rds objects
#   - CSV summaries for supplement/QA
# ============================================================

# ---- EDIT THESE PATHS ----
in_dir  <- "/part1/02_pathologist_annot"   # output from 1.1
out_dir <- "/part1/03_qc_spotsweeper"      # output of 1.2
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- GLOBAL CUT-OFFS (your current defaults) ----
nCount_cutoff   <- 200
nFeature_cutoff <- 200

# -----------------------------
# Helpers
# -----------------------------
extract_core_10x_barcode <- function(seurat_barcodes) {
  out <- sub(".*_(.*-1)$", "\\1", seurat_barcodes)
  out
}

seurat_to_spe <- function(seu) {
  coords <- GetTissueCoordinates(seu, scale = "hires")[, c("x", "y")]
  SpatialExperiment(
    assays = list(counts = seu@assays$Spatial$counts),
    colData = seu@meta.data,
    spatialCoords = as.matrix(coords)
  )
}

transfer_spe_coldata_to_seurat <- function(seu, spe) {
  spe_meta <- as.data.frame(colData(spe))
  # ensure same order by cell names
  if (!all(rownames(spe_meta) == colnames(seu))) {
    spe_meta <- spe_meta[colnames(seu), , drop = FALSE]
  }
  for (col in colnames(spe_meta)) {
    seu@meta.data[[col]] <- spe_meta[[col]]
  }
  seu
}

# -----------------------------
# Load all inputs automatically
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in: ", in_dir)

sample_ids <- sub("\\.rds$|\\.RDS$", "", basename(rds_files))
names(rds_files) <- sample_ids

message("Found ", length(rds_files), " samples")

# -----------------------------
# Run SpotSweeper QC per sample
# -----------------------------
outlier_summary <- data.frame(
  Sample = sample_ids,
  Total_Spots = NA_integer_,
  Outlier_Spots = NA_integer_,
  Percent_Outliers = NA_real_,
  Spots_Kept_After_Sweeper = NA_integer_,
  stringsAsFactors = FALSE
)

filtered_summary <- data.frame(
  Sample = sample_ids,
  Original_Spots = NA_integer_,
  After_SpotSweeper = NA_integer_,
  After_GlobalCutoff = NA_integer_,
  Percent_Retained = NA_real_,
  stringsAsFactors = FALSE
)

for (sid in sample_ids) {
  message("\n--- ", sid, " ---")
  seu <- readRDS(rds_files[[sid]])
  
  # convert to SpatialExperiment
  spe <- seurat_to_spe(seu)
  
  # mito genes
  gene_names <- rownames(assay(spe, "counts"))
  mito <- gene_names[grepl("^MT-", gene_names)]
  
  spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = mito))
  
  # SpotSweeper local outliers
  spe <- SpotSweeper::localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
  spe <- SpotSweeper::localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
  spe <- SpotSweeper::localOutliers(spe, metric = "subsets_Mito_percent", direction = "higher", log = FALSE)
  
  spe$local_outliers <- as.logical(spe$sum_outliers) |
    as.logical(spe$detected_outliers) |
    as.logical(spe$subsets_Mito_percent_outliers)
  
  n_total <- ncol(seu)
  n_out   <- sum(spe$local_outliers %in% TRUE, na.rm = TRUE)
  
  outlier_summary[outlier_summary$Sample == sid, "Total_Spots"] <- n_total
  outlier_summary[outlier_summary$Sample == sid, "Outlier_Spots"] <- n_out
  outlier_summary[outlier_summary$Sample == sid, "Percent_Outliers"] <- round(100 * n_out / n_total, 2)
  outlier_summary[outlier_summary$Sample == sid, "Spots_Kept_After_Sweeper"] <- n_total - n_out
  
  # transfer QC flags to Seurat metadata
  seu <- transfer_spe_coldata_to_seurat(seu, spe)
  
  # filter: remove local outliers
  keep1 <- is.na(seu$local_outliers) | seu$local_outliers == FALSE
  seu1 <- subset(seu, cells = colnames(seu)[keep1])
  
  # filter: global cutoffs
  keep2 <- (seu1$nCount_Spatial > nCount_cutoff) & (seu1$nFeature_Spatial > nFeature_cutoff)
  seu2 <- subset(seu1, cells = colnames(seu1)[keep2])
  
  # fill filtered summary
  filtered_summary[filtered_summary$Sample == sid, "Original_Spots"] <- n_total
  filtered_summary[filtered_summary$Sample == sid, "After_SpotSweeper"] <- ncol(seu1)
  filtered_summary[filtered_summary$Sample == sid, "After_GlobalCutoff"] <- ncol(seu2)
  filtered_summary[filtered_summary$Sample == sid, "Percent_Retained"] <- round(100 * ncol(seu2) / n_total, 1)
  
  message("Kept ", ncol(seu2), " / ", n_total, " spots after all filtering.")
  
  saveRDS(seu2, file.path(out_dir, paste0(sid, ".rds")))
}

# write summaries
write.csv(outlier_summary, file.path(out_dir, "SpotSweeper_outlier_summary.csv"), row.names = FALSE)
write.csv(filtered_summary, file.path(out_dir, "Supplementary_SpotFiltering_Summary.csv"), row.names = FALSE)

message("\nDone.")
message("Saved filtered objects to: ", out_dir)
message("Saved: SpotSweeper_outlier_summary.csv")
message("Saved: Supplementary_SpotFiltering_Summary.csv")
