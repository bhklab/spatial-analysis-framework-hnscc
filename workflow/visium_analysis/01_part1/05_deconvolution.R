suppressPackageStartupMessages({
  library(Seurat)
  library(CARD)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(scuttle)
})

# ============================================================
# Part 1.4 — Deconvolution (CARD)
# Input:  Seurat objects after QC (from Part 1.2)
# Output: CARD objects + proportion tables
# ============================================================

# ---- EDIT THESE PATHS ----
in_dir <- "/part1/03_qc_spotsweeper"           # output from 1.2
out_dir <- "/part1/04_deconvolution"          # new output folder
ref_dir <- "/part1/04_sc_reference/refs_diet" # from Part 1.3 (diet refs)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load diet references (from your 1.3 outputs) ----
ref_paths <- list(
  choi  = file.path(ref_dir, "choi_ref_diet.rds"),
  quah  = file.path(ref_dir, "quah_ref_diet.rds"),
  puram = file.path(ref_dir, "puram_ref_diet.rds")
)

# -----------------------------
# Helpers
# -----------------------------
seurat_to_card_inputs <- function(seu) {
  # CARD needs: spatial_count, spatial_location
  counts <- GetAssayData(seu, assay = "Spatial", slot = "counts")
  coords <- GetTissueCoordinates(seu, scale = "hires")[, c("x", "y")]
  coords <- as.data.frame(coords)
  coords$spot <- rownames(coords)
  
  # Ensure same ordering
  coords <- coords[colnames(seu), , drop = FALSE]
  rownames(coords) <- coords$spot
  
  list(spatial_count = counts, spatial_location = coords)
}

prep_reference_for_card <- function(ref_seu, ct_col = "cell_type", sample_col = "orig.ident") {
  stopifnot(ct_col %in% colnames(ref_seu@meta.data))
  stopifnot(sample_col %in% colnames(ref_seu@meta.data))
  
  sc_count <- GetAssayData(ref_seu, assay = "RNA", slot = "counts")
  sc_meta <- ref_seu@meta.data
  sc_meta$cell_type <- sc_meta[[ct_col]]
  sc_meta$orig.ident <- sc_meta[[sample_col]]
  
  list(sc_count = sc_count, sc_meta = sc_meta)
}

run_card_for_one_sample <- function(spatial_count, spatial_location, sc_count, sc_meta,
                                    minCountGene = 100, minCountSpot = 5) {
  card_obj <- createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "cell_type",
    ct.select = unique(sc_meta$cell_type),
    sample.varname = "orig.ident",
    minCountGene = minCountGene,
    minCountSpot = minCountSpot
  )
  CARD_deconvolution(CARD_object = card_obj)
}

# -----------------------------
# Main
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in: ", in_dir)

sample_ids <- sub("\\.rds$|\\.RDS$", "", basename(rds_files))
names(rds_files) <- sample_ids

for (ref_name in names(ref_paths)) {
  message("\n==============================")
  message("Reference: ", ref_name)
  message("==============================")
  
  ref_path <- ref_paths[[ref_name]]
  if (!file.exists(ref_path)) stop("Missing reference: ", ref_path)
  
  ref_seu <- readRDS(ref_path)
  
  # OPTIONAL: apply your reference filtering here (e.g., remove Stress/Cycling)
  # Idents(ref_seu) <- "cell_type"
  # ref_seu <- subset(ref_seu, idents = c("Stress cells", "Cycling cells"), invert = TRUE)
  
  ref <- prep_reference_for_card(ref_seu, ct_col = "cell_type", sample_col = "orig.ident")
  
  ref_out_dir <- file.path(out_dir, ref_name)
  dir.create(ref_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sid in sample_ids) {
    message("\n--- ", sid, " (", ref_name, ") ---")
    seu <- readRDS(rds_files[[sid]])
    
    inp <- seurat_to_card_inputs(seu)
    
    card_res <- run_card_for_one_sample(
      spatial_count = inp$spatial_count,
      spatial_location = inp$spatial_location,
      sc_count = ref$sc_count,
      sc_meta = ref$sc_meta
    )
    
    saveRDS(card_res, file.path(ref_out_dir, paste0(sid, "_CARD.rds")))
    
    # also save proportions if present
    if (!is.null(card_res@Proportion_CELL_TYPE)) {
      write.csv(card_res@Proportion_CELL_TYPE,
                file.path(ref_out_dir, paste0(sid, "_CARD_proportions.csv")))
    }
  }
}

message("\nDone. Outputs in: ", out_dir)
