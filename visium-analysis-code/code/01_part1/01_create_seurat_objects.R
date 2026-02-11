suppressPackageStartupMessages({
  library(Seurat)
  library(DropletUtils)
  library(yaml)
})


# -----------------------------
# Helpers
# -----------------------------
ensure_h5_from_filtered_matrix <- function(sample_dir,
                                           matrix_subdir = "filtered_feature_bc_matrix",
                                           h5_name = "filtered_feature_bc_matrix.h5") {
  matrix_dir <- file.path(sample_dir, matrix_subdir)
  h5_path    <- file.path(sample_dir, h5_name)
  
  if (file.exists(h5_path)) {
    message("  - found existing: ", h5_path)
    return(invisible(h5_path))
  }
  
  if (!dir.exists(matrix_dir)) {
    stop("Missing matrix directory: ", matrix_dir)
  }
  
  message("  - building h5 from matrix dir: ", matrix_dir)
  mat <- Read10X(matrix_dir)
  if (is.list(mat)) mat <- mat[[1]]
  
  DropletUtils::write10xCounts(
    h5_path,
    mat,
    type = "HDF5",
    version = "3",
    overwrite = TRUE,
    gene.id = rownames(mat),
    gene.symbol = rownames(mat)
  )
  
  if (!file.exists(h5_path)) stop("Failed to create: ", h5_path)
  invisible(h5_path)
}

# -----------------------------
# Main
# -----------------------------
cfg <- yaml::read_yaml("configs/config.yml")
samples <- read.csv("configs/samples_part1.csv", stringsAsFactors = FALSE)

spaceranger_root <- cfg$paths$spaceranger_root
out_dir <- file.path(cfg$paths$output_root, "part1/01_seurat_raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assay <- cfg$run$assay %||% "Spatial"
slice <- cfg$run$slice %||% "spatial"
filter_matrix <- cfg$run$filter_matrix %||% TRUE

# tiny infix helper like in rlang, avoids importing more packages
`%||%` <- function(x, y) if (is.null(x)) y else x

for (i in seq_len(nrow(samples))) {
  sid  <- samples$sample_id[i]
  sdir <- samples$spaceranger_dir[i]
  tp   <- if ("timepoint" %in% names(samples)) samples$timepoint[i] else NA
  
  message("\n[", i, "/", nrow(samples), "] ", sid, "  (", sdir, ")")
  
  sample_dir <- file.path(spaceranger_root, sdir)
  if (!dir.exists(sample_dir)) stop("Missing sample_dir: ", sample_dir)
  
  # ensure filtered_feature_bc_matrix.h5 exists (your current workflow)
  ensure_h5_from_filtered_matrix(sample_dir)
  
  # Build Seurat object
  obj <- Seurat::Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = assay,
    slice = slice,
    filter.matrix = filter_matrix,
    image = NULL
  )
  
  # Add consistent metadata
  obj$sample_id <- sid
  obj$timepoint <- tp
  obj$spaceranger_dir <- sdir
  obj$orig.ident <- sid
  
  saveRDS(obj, file.path(out_dir, paste0(sid, ".rds")))
  message("  - saved: ", file.path(out_dir, paste0(sid, ".rds")))
}

message("\nDone. Outputs in: ", out_dir)
