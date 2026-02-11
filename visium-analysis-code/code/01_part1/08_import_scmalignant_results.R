suppressPackageStartupMessages({
  library(Seurat)
})

# ============================================================
# Part 1.5 — Import scMalignantFinder results into Seurat objects
# Reads:
#   - <sample>_malignant_results.csv
# Adds columns to meta.data:
#   - malignancy_probability
#   - scMalignantFinder_prediction
# ============================================================

# ---- EDIT PATHS ----
seurat_dir  <- "/part1/04_deconvolution_or_latest_rds"  # same objects used for export
results_dir <- "/part1/06_scmf_output"                   # where python wrote *_malignant_results.csv
out_dir     <- "/part1/07_malignant_annot"               # output annotated objects

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(seurat_dir, pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in: ", seurat_dir)

for (f in rds_files) {
  sample <- sub("\\.rds$|\\.RDS$", "", basename(f))
  message("Annotating: ", sample)
  
  seu <- readRDS(f)
  
  csv_file <- file.path(results_dir, paste0(sample, "_malignant_results.csv"))
  if (!file.exists(csv_file)) {
    message("  [WARN] Missing results for ", sample, ": ", csv_file)
    next
  }
  
  # rows are barcodes; columns include malignancy_probability + prediction
  res <- read.csv(csv_file, row.names = 1, stringsAsFactors = FALSE)
  
  # enforce matching by barcode
  common <- intersect(rownames(res), colnames(seu))
  if (length(common) == 0) {
    message("  [WARN] No barcode overlap for ", sample, ". Skipping.")
    next
  }
  
  res <- res[common, , drop = FALSE]
  seu <- AddMetaData(seu, metadata = res)
  
  saveRDS(seu, file.path(out_dir, paste0(sample, "_Annotated.rds")))
  message("  saved: ", file.path(out_dir, paste0(sample, "_Annotated.rds")))
}

message("Done. Annotated Seurat objects in: ", out_dir)
