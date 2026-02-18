suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

# ============================================================
# Part 1.3 — Build single-cell reference objects (Choi / Quah / Puram)
# Outputs are saved for downstream deconvolution (Part 1.5).
# ============================================================

# ---- EDIT THESE PATHS ----
in_choi_counts <- "/GSE181919_UMI_counts.txt"
in_choi_meta   <- "/GSE181919_Barcode_metadata.txt"

in_quah_rds    <- "/Quah_reference.rds"
in_puram_rds   <- "/Puram_reference.rds"

out_dir <- "/part1/04_sc_reference"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
safe_diet_seurat <- function(obj) {
  # DietSeurat sometimes drops needed things; we keep counts + meta
  obj@assays$RNA@data <- NULL
  obj@assays$RNA@scale.data <- NULL
  obj
}

standard_preprocess <- function(obj, dims = 1:30, do_integrate = FALSE, batch_col = NULL) {
  # Minimal, reproducible reference preprocessing
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  
  if (do_integrate) {
    stopifnot(!is.null(batch_col), batch_col %in% colnames(obj@meta.data))
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[batch_col]][, 1])
    obj <- SCTransform(obj, verbose = FALSE)
    obj <- SelectIntegrationFeatures(obj, nfeatures = 3000)
    obj <- PrepSCTIntegration(obj)
    anchors <- FindIntegrationAnchors(obj, normalization.method = "SCT")
    obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    DefaultAssay(obj) <- "integrated"
    obj <- RunPCA(obj)
    obj <- RunUMAP(obj, dims = dims)
  } else {
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = dims)
  }
  
  obj
}

# -----------------------------
# 1) Choi (GSE181919) from text matrices
# -----------------------------
message("Building Choi reference...")
choi_counts <- read.delim2(in_choi_counts, check.names = FALSE)
choi_meta   <- read.delim2(in_choi_meta, check.names = FALSE)

# your original code replaces "." -> "-" to match metadata rownames
colnames(choi_counts) <- sub("\\.", "-", colnames(choi_counts))

# ensure metadata rownames are barcodes
if (is.null(rownames(choi_meta)) || any(rownames(choi_meta) == "")) {
  # if barcode column exists, use it (change "barcode" if needed)
  if ("barcode" %in% colnames(choi_meta)) rownames(choi_meta) <- choi_meta$barcode
}

# sanity check alignment
stopifnot(all(colnames(choi_counts) %in% rownames(choi_meta)))

choi <- CreateSeuratObject(counts = choi_counts, meta.data = choi_meta)

# keep LN + primary if your meta has tissue.type with CA/LN
if ("tissue.type" %in% colnames(choi@meta.data)) {
  choi <- subset(choi, subset = tissue.type %in% c("CA", "LN"))
}

# If you were integrating by patient.id (your code splits by `patient.id` / `pa`)
batch_col <- NULL
if ("patient.id" %in% colnames(choi@meta.data)) batch_col <- "patient.id"
if ("pa" %in% colnames(choi@meta.data)) batch_col <- "pa"

if (!is.null(batch_col)) {
  choi <- standard_preprocess(choi, dims = 1:30, do_integrate = TRUE, batch_col = batch_col)
} else {
  choi <- standard_preprocess(choi, dims = 1:30, do_integrate = FALSE)
}

saveRDS(choi, file.path(out_dir, "choi_ref.rds"))

# -----------------------------
# 2) Quah reference (already RDS)
# -----------------------------
message("Building Quah reference...")
quah <- readRDS(in_quah_rds)

# integrate if it has a batch column; otherwise do minimal preprocess
batch_col <- NULL
if ("patientID" %in% colnames(quah@meta.data)) batch_col <- "patientID"
if (!is.null(batch_col)) {
  quah <- standard_preprocess(quah, dims = 1:30, do_integrate = TRUE, batch_col = batch_col)
} else {
  quah <- standard_preprocess(quah, dims = 1:30, do_integrate = FALSE)
}

# your original script removes “Salivary Gland”
if ("cell_type" %in% colnames(quah@meta.data)) {
  quah <- subset(quah, subset = cell_type != "Salivary Gland")
}

saveRDS(quah, file.path(out_dir, "quah_ref.rds"))

# -----------------------------
# 3) Puram reference (already RDS)
# -----------------------------
message("Building Puram reference...")
puram <- readRDS(in_puram_rds)

batch_col <- NULL
if ("patient_id" %in% colnames(puram@meta.data)) batch_col <- "patient_id"
if (!is.null(batch_col)) {
  puram <- standard_preprocess(puram, dims = 1:30, do_integrate = TRUE, batch_col = batch_col)
} else {
  puram <- standard_preprocess(puram, dims = 1:30, do_integrate = FALSE)
}

saveRDS(puram, file.path(out_dir, "puram_ref.rds"))

# -----------------------------
# 4) Diet versions (for memory-efficient deconvolution)
# -----------------------------
message("Saving diet references...")
diet_dir <- file.path(out_dir, "refs_diet")
dir.create(diet_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(safe_diet_seurat(choi),  file.path(diet_dir, "choi_ref_diet.rds"))
saveRDS(safe_diet_seurat(quah),  file.path(diet_dir, "quah_ref_diet.rds"))
saveRDS(safe_diet_seurat(puram), file.path(diet_dir, "puram_ref_diet.rds"))

message("Done. References in: ", out_dir)
