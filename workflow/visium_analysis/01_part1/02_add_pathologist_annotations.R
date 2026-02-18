suppressPackageStartupMessages({
  library(Seurat)
})

# ============================================================
# Part 1.1 — Add Pathologist Annotations to Visium Seurat Objects
# Inputs:
#   - in_dir:  folder containing *.rds Seurat objects (from Part 1.0)
#   - anno_dir: folder containing per-sample CSV annotations (e.g., S1_Pre.csv)
# Output:
#   - out_dir: folder containing annotated *.rds objects
#
# Required columns in annotation CSV:
#   - Barcode
#   - Patho_Anno_Jinsu  
# ============================================================

# ---- EDIT THESE THREE PATHS ----
in_dir  <- "/part1/01_seurat_raw"                 # output from 1.0
anno_dir <- "/pathologist_annotation_csv"        # has <sample_id>.csv
out_dir <- "/part1/02_pathologist_annot"          # output of 1.1

# ---- SETTINGS ----
anno_barcode_col <- "Barcode"
anno_value_col   <- "Patho_Anno_Jinsu"     # column in CSV
meta_col_name    <- "Pathologist_Annotation"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: robust barcode alignment
# -----------------------------
standardize_ann_barcodes <- function(x) {
  x <- trimws(x)
  # ensure "-1" suffix
  ifelse(grepl("-1$", x), x, paste0(x, "-1"))
}

extract_core_10x_barcode <- function(seurat_barcodes) {
  # keep core 10x barcode at the end
  sub(".*_(.*-1)$", "\\1", seurat_barcodes)
}

add_pathologist_annotation <- function(seu_obj, anno_df,
                                       anno_barcode_col, anno_value_col,
                                       meta_col_name) {
  if (!all(c(anno_barcode_col, anno_value_col) %in% colnames(anno_df))) {
    stop("Missing required columns: ", paste(setdiff(c(anno_barcode_col, anno_value_col), colnames(anno_df)), collapse=", "))
  }
  
  seurat_barcodes <- colnames(seu_obj)
  core_barcodes   <- extract_core_10x_barcode(seurat_barcodes)
  
  anno_barcodes <- standardize_ann_barcodes(anno_df[[anno_barcode_col]])
  anno_values   <- trimws(as.character(anno_df[[anno_value_col]]))
  
  # map: core barcode -> full seurat barcode
  barcode_map <- setNames(seurat_barcodes, core_barcodes)
  # map: anno barcode -> annotation value
  anno_vec <- setNames(anno_values, anno_barcodes)
  
  matched <- intersect(names(anno_vec), names(barcode_map))
  if (length(matched) == 0) {
    return(list(seu = seu_obj, n_matched = 0L))
  }
  
  # annotation values named by full seurat barcode
  annotated_values <- setNames(anno_vec[matched], barcode_map[matched])
  
  seu_obj <- AddMetaData(
    object = seu_obj,
    metadata = annotated_values,
    col.name = meta_col_name
  )
  
  list(seu = seu_obj, n_matched = length(matched))
}

# -----------------------------
# Main: iterate over all RDS files
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in: ", in_dir)

message("Found ", length(rds_files), " input objects in: ", in_dir)

summary_df <- data.frame(
  sample_id = character(),
  n_spots = integer(),
  n_matched = integer(),
  anno_file_found = logical(),
  stringsAsFactors = FALSE
)

for (f in rds_files) {
  sample_id <- sub("\\.rds$|\\.RDS$", "", basename(f))
  message("\n--- ", sample_id, " ---")
  
  anno_file <- file.path(anno_dir, paste0(sample_id, ".csv"))
  if (!file.exists(anno_file)) {
    message("❌ Missing annotation CSV: ", anno_file)
    seu_obj <- readRDS(f)
    summary_df <- rbind(summary_df, data.frame(
      sample_id = sample_id,
      n_spots = ncol(seu_obj),
      n_matched = 0L,
      anno_file_found = FALSE
    ))
    next
  }
  
  seu_obj <- readRDS(f)
  anno_df <- read.csv(anno_file, stringsAsFactors = FALSE)
  
  res <- add_pathologist_annotation(
    seu_obj, anno_df,
    anno_barcode_col = anno_barcode_col,
    anno_value_col   = anno_value_col,
    meta_col_name    = meta_col_name
  )
  
  message("✅ Matched barcodes: ", res$n_matched, " / ", ncol(seu_obj))
  
  saveRDS(res$seu, file.path(out_dir, paste0(sample_id, ".rds")))
  
  summary_df <- rbind(summary_df, data.frame(
    sample_id = sample_id,
    n_spots = ncol(seu_obj),
    n_matched = res$n_matched,
    anno_file_found = TRUE
  ))
}

# Save a small run summary (very publication-friendly)
write.csv(summary_df, file.path(out_dir, "pathologist_annotation_summary.csv"), row.names = FALSE)

message("\nDone. Annotated objects saved to: ", out_dir)
message("Summary written to: ", file.path(out_dir, "pathologist_annotation_summary.csv"))
