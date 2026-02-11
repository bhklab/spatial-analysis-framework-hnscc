suppressPackageStartupMessages({
  library(Seurat)
})

# ============================================================
# Part 1.5 — Export Visium Seurat objects for scMalignantFinder
# Exports per-sample:
#   - <sample>_counts.csv  (genes x spots)
#   - <sample>_metadata.csv (spots x metadata)
# ============================================================

# ---- EDIT PATHS ----
in_dir  <- "/part1/04_deconvolution_or_latest_rds"  # folder containing *.rds Seurat objects
out_dir <- "/part1/05_scmf_input"                  # where to write scMF inputs

assay_name <- "Spatial"   # usually Spatial for Visium

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(in_dir, pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in: ", in_dir)

for (f in rds_files) {
  sample_name <- sub("\\.rds$|\\.RDS$", "", basename(f))
  message("Exporting: ", sample_name)
  
  seu <- readRDS(f)
  
  if (!assay_name %in% names(seu@assays)) {
    stop("Assay '", assay_name, "' not found in ", sample_name,
         ". Available assays: ", paste(names(seu@assays), collapse = ", "))
  }
  
  # counts: genes x spots
  counts <- GetAssayData(seu, assay = assay_name, slot = "counts")
  counts_mat <- as.matrix(counts)
  
  # metadata: spots x meta
  meta <- seu@meta.data
  
  # ensure matching order
  meta <- meta[colnames(seu), , drop = FALSE]
  
  write.csv(counts_mat, file = file.path(out_dir, paste0(sample_name, "_counts.csv")))
  write.csv(meta,       file = file.path(out_dir, paste0(sample_name, "_metadata.csv")))
  
  message("  - wrote: ", sample_name, "_counts.csv + _metadata.csv")
}

message("Done. scMF inputs in: ", out_dir)
