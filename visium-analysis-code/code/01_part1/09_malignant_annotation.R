suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Part 1.6 — Final malignant annotation + broad celltype label
# Input:
#   - in_dir: per-sample Seurat objects AFTER scMalignantFinder import
#   - sample_info.csv + patient_info.csv for clinical metadata
# Output:
#   - out_dir/per_sample/<sample>.rds
#   - out_dir/AllST.rds
# ============================================================

# ---- EDIT THESE PATHS ----
in_dir <- "/part1/07_malignant_annot"      # your per-sample *_Annotated.rds folder
sample_info_csv  <- "/sample_info.csv"
patient_info_csv <- "/patient_info.csv"
out_dir <- "/part1/08_final_malignant_annotation"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "per_sample"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
make_sample_prefix <- function(id_patient, id_sample) {
  # matches your logic:
  # LIB_04 -> LIB04 ; INS_12 -> S12 ; else take prefix from id_sample
  if (grepl("^LIB", id_patient)) {
    paste0("LIB", str_pad(as.numeric(str_extract(id_patient, "\\d+$")), width = 2, pad = "0"))
  } else if (grepl("^INS", id_patient)) {
    paste0("S", as.numeric(str_extract(id_patient, "\\d+$")))
  } else {
    sub("_.*", "", id_sample)
  }
}

make_timepoint <- function(sample_timing) {
  ifelse(grepl("BL|Pre", sample_timing, ignore.case = TRUE), "Pre", "Post")
}

# Assign broad “Other” label based on deconv proportion columns (CARD/RCTD)
assign_other_dominant <- function(meta, final_col = "Final_Malignant_Annotation") {
  # find all deconvolution columns except the “Dominant” categorical ones
  all_deconv_cols <- grep("^(CARD|RCTD)_", colnames(meta), value = TRUE)
  dominant_cols <- grep("Dominant", all_deconv_cols, ignore.case = TRUE, value = TRUE)
  
  # keep numeric proportion columns only
  candidate <- setdiff(all_deconv_cols, dominant_cols)
  candidate <- candidate[sapply(meta[, candidate, drop = FALSE], is.numeric)]
  
  if (length(candidate) == 0) {
    meta$Other_Dominant_Celltype <- NA_character_
    return(meta)
  }
  
  # Standardize names from column headers (your original approach)
  raw_labels <- sub("^(CARD|RCTD)_(Puram|Quah|Choi)[._]", "", candidate)
  raw_labels <- gsub("[._]cells$", "", raw_labels)
  
  standardized <- tolower(raw_labels)
  standardized <- gsub("^(b|t|dendritic|endothelial|fibroblast|macrophage).*", "\\1", standardized)
  
  names(standardized) <- candidate
  
  other_idx <- which(meta[[final_col]] != "Malignant" | is.na(meta[[final_col]]))
  meta$Other_Dominant_Celltype <- NA_character_
  
  if (length(other_idx) > 0) {
    assign_one <- function(row) {
      scores <- as.numeric(row[candidate])
      names(scores) <- standardized
      scores_agg <- tapply(scores, names(scores), mean, na.rm = TRUE)
      names(which.max(scores_agg))
    }
    meta$Other_Dominant_Celltype[other_idx] <- apply(meta[other_idx, candidate, drop = FALSE], 1, assign_one)
  }
  
  pretty <- c(
    "b" = "B cells",
    "t" = "T cells",
    "dendritic" = "Dendritic cells",
    "endothelial" = "Endothelial cells",
    "fibroblast" = "Fibroblasts",
    "macrophage" = "Macrophages"
  )
  meta$Other_Dominant_Celltype <- unname(pretty[meta$Other_Dominant_Celltype])
  
  meta
}

# -----------------------------
# 1) Load all per-sample objects automatically
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in: ", in_dir)

sample_ids <- sub("\\.rds$|\\.RDS$", "", basename(rds_files))
names(rds_files) <- sample_ids

message("Found ", length(rds_files), " samples in: ", in_dir)

seurat_objects <- lapply(rds_files, readRDS)

# -----------------------------
# 2) Clinical metadata merge (sample_info + patient_info)
# -----------------------------
clinical_data <- read.csv(sample_info_csv, stringsAsFactors = FALSE)
patient_data  <- read.csv(patient_info_csv, stringsAsFactors = FALSE)

clinical_data$Sample_Prefix <- mapply(make_sample_prefix, clinical_data$id_patient, clinical_data$id_sample)
clinical_data$Timepoint <- make_timepoint(clinical_data$sample_timing)
clinical_data$Seurat_Name <- paste0(clinical_data$Sample_Prefix, "_", clinical_data$Timepoint)

patient_data$Sample_Prefix <- ifelse(
  grepl("^LIB_", patient_data$id_study),
  paste0("LIB", str_pad(as.numeric(str_extract(patient_data$id_study, "\\d+$")), width = 2, pad = "0")),
  patient_data$id_study
)

clinical_merged <- dplyr::left_join(clinical_data, patient_data, by = "Sample_Prefix")

# attach metadata to each object
for (sid in names(seurat_objects)) {
  obj <- seurat_objects[[sid]]
  row <- clinical_merged[clinical_merged$Seurat_Name == sid, , drop = FALSE]
  
  if (nrow(row) == 1) {
    cols_to_add <- setdiff(colnames(row),
                           c("Seurat_Name", "Sample_Prefix", "Timepoint", "id_sample", "id_patient", "id_study"))
    for (col in cols_to_add) obj[[col]] <- row[[col]]
  } else {
    warning("No (or multiple) clinical rows matched for: ", sid)
  }
  
  seurat_objects[[sid]] <- obj
}

# -----------------------------
# 3) Rename cells + merge into AllST
# -----------------------------
for (sid in names(seurat_objects)) {
  seurat_objects[[sid]] <- RenameCells(seurat_objects[[sid]], add.cell.id = sid)
  seurat_objects[[sid]]$V <- Cells(seurat_objects[[sid]])
}

options(future.globals.maxSize = 30 * 1024^3)
AllST <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])

# -----------------------------
# 4) Final malignant rule (your exact logic)
# -----------------------------
stopifnot("Pathologist_Annotation" %in% colnames(AllST@meta.data))
stopifnot("scMalignantFinder_prediction" %in% colnames(AllST@meta.data))

deconv_cols <- c(
  "CARD_Choi_Dominant_celltype",
  "CARD_Puram_Dominant_celltype",
  "CARD_Quah_Dominant_celltype",
  "RCTD_Puram_Dominant_Celltype",
  "RCTD_Choi_Dominant_Celltype",
  "RCTD_Quah_Dominant_Celltype"
)
present_deconv_cols <- deconv_cols[deconv_cols %in% colnames(AllST@meta.data)]
if (length(present_deconv_cols) == 0) {
  warning("No deconv dominant columns found; deconv consensus will be skipped.")
  deconv_malignant_count <- rep(0L, ncol(AllST))
} else {
  deconv_malignant_count <- rowSums(AllST@meta.data[, present_deconv_cols, drop = FALSE] == "malignant_cells", na.rm = TRUE)
}

AllST$Final_Malignant_Annotation <- ifelse(
  AllST$Pathologist_Annotation == "SCC" &
    (AllST$scMalignantFinder_prediction == "Malignant") &
    (deconv_malignant_count >= 4),
  "Malignant", "Other"
)

# -----------------------------
# 5) Broad_CT for non-malignant spots (from numeric deconv proportions)
# -----------------------------
AllST@meta.data <- assign_other_dominant(AllST@meta.data, final_col = "Final_Malignant_Annotation")

AllST$Broad_CT <- ifelse(
  AllST$Final_Malignant_Annotation == "Malignant",
  "Malignant cells",
  AllST$Other_Dominant_Celltype
)

# -----------------------------
# 6) Split back to per-sample and save
# -----------------------------
AllST$ID <- sub("_[^_]+$", "", colnames(AllST))  # same logic as your original
samples_split <- SplitObject(AllST, split.by = "ID")

for (sid in names(samples_split)) {
  saveRDS(samples_split[[sid]], file.path(out_dir, "per_sample", paste0(sid, ".rds")))
}
saveRDS(AllST, file.path(out_dir, "AllST.rds"))

message("Done.")
message("Saved AllST: ", file.path(out_dir, "AllST.rds"))
message("Saved per-sample objects in: ", file.path(out_dir, "per_sample"))
