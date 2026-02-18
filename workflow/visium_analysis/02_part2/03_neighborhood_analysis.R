suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(spdep)
})

# ============================================================
# Part 2.2a — Neighborhood analysis (kNN=6)
#
# Input:
#   RM_ST_with_nonmalignant_refined_labels.rds  (from Part 2.1 refactor)
#   Must contain a label column used for neighborhood counting.
#
# Output:
#   combined_neighborhood_counts.csv
#   Sample_Summary_Table.csv
#
# ============================================================

# ---- EDIT PATHS ----
in_rds  <- "/part2/02_nonmalignant_characterization/RM_ST_with_nonmalignant_refined_labels.rds"
out_dir <- "/part2/03_neighborhood_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_counts_csv <- file.path(out_dir, "combined_neighborhood_counts.csv")
out_summary_csv <- file.path(out_dir, "Sample_Summary_Table.csv")

# ---- PARAMETERS ----
label_col <- "Spatial_CT_refined"   # created in Part 2.1 refactor
k_neighbors <- 6

# -----------------------------
# Spatial coordinate extraction
# -----------------------------
extract_spatial_coords <- function(seurat_obj) {
  # Choose the first image that exists
  if (length(seurat_obj@images) == 0) stop("No spatial images found in object.")
  img_name <- names(seurat_obj@images)[[1]]
  
  # Try VisiumV2 style
  coords <- tryCatch({
    c <- GetTissueCoordinates(seurat_obj, image = img_name)
    if (all(c("x","y") %in% colnames(c))) {
      c <- c[, c("x","y"), drop = FALSE]
      colnames(c) <- c("col","row")
    }
    c
  }, error = function(e) NULL)
  
  # Fall back to @coordinates (VisiumV1-style)
  if (is.null(coords)) {
    img <- seurat_obj@images[[img_name]]
    if (!("coordinates" %in% slotNames(img))) stop("No coordinates slot in spatial image.")
    coords <- img@coordinates
    
    # common Visium V1 cases
    if (all(c("tissue","row","col") %in% colnames(coords))) {
      coords <- coords[coords$tissue == 1, c("col","row"), drop = FALSE]
    } else if (all(c("array_row","array_col") %in% colnames(coords))) {
      coords <- coords[, c("array_col","array_row"), drop = FALSE]
      colnames(coords) <- c("col","row")
    } else if (all(c("imagecol","imagerow") %in% colnames(coords))) {
      coords <- coords[, c("imagecol","imagerow"), drop = FALSE]
      colnames(coords) <- c("col","row")
    } else {
      stop("Unrecognized coordinate columns in image@coordinates.")
    }
  }
  
  # Ensure rownames are spot IDs
  if (is.null(rownames(coords))) stop("Coordinates have no rownames (spot IDs).")
  coords
}

# -----------------------------
# Neighbor computation (kNN=6)
# -----------------------------
compute_neighbors_knn <- function(seurat_obj, annotation_col, k = 6) {
  coords <- extract_spatial_coords(seurat_obj)
  
  # Keep only spots that exist in the object
  keep <- intersect(rownames(coords), colnames(seurat_obj))
  coords <- coords[keep, , drop = FALSE]
  
  coords_mat <- as.matrix(coords[, c("col","row"), drop = FALSE])
  
  # kNN neighborhood
  nn <- spdep::knn2nb(spdep::knearneigh(coords_mat, k = k))
  
  # Store neighbor barcodes + neighbor labels
  md <- seurat_obj@meta.data
  
  md[keep, "neighbor_bcs"] <- vapply(seq_along(nn), function(i) {
    paste(rownames(coords_mat)[nn[[i]]], collapse = ",")
  }, character(1))
  
  md[keep, "neighbor_clust_anno"] <- vapply(seq_along(nn), function(i) {
    nbs <- rownames(coords_mat)[nn[[i]]]
    paste(md[nbs, annotation_col], collapse = ",")
  }, character(1))
  
  seurat_obj@meta.data <- md
  seurat_obj
}

# -----------------------------
# Count malignant neighbors
# -----------------------------
count_neighbors <- function(obj, label_col) {
  stopifnot(label_col %in% colnames(obj@meta.data))
  stopifnot(all(c("neighbor_clust_anno","neighbor_bcs") %in% colnames(obj@meta.data)))
  
  md <- obj@meta.data[, c(label_col, "neighbor_clust_anno", "neighbor_bcs"), drop = FALSE]
  colnames(md)[1] <- "primary"
  
  # explode comma-separated neighbor label list
  neigh_df <- data.frame(
    primary = rep(md$primary, lengths(strsplit(md$neighbor_clust_anno, ","))),
    neighbor = unlist(strsplit(md$neighbor_clust_anno, ",")),
    neighbor_barcodes = unlist(strsplit(md$neighbor_bcs, ","))
  )
  
  # malignant set = Cancer_C*
  malignant_levels <- unique(md$primary)
  malignant_set <- malignant_levels[grepl("^Cancer_C", malignant_levels)]
  
  # around malignant, exclude malignant neighbors
  neigh_df %>%
    filter(primary %in% malignant_set) %>%
    filter(!neighbor %in% malignant_set) %>%
    group_by(primary, neighbor) %>%
    summarise(`#neighbors` = n(), .groups = "drop")
}

# ============================================================
# Main
# ============================================================
RM_ST <- readRDS(in_rds)
stopifnot(label_col %in% colnames(RM_ST@meta.data))
RM_ST <- tryCatch(JoinLayers(RM_ST), error = function(e) RM_ST)

# Split by sample ID (must exist)
stopifnot("ID" %in% colnames(RM_ST@meta.data))
samples <- SplitObject(RM_ST, split.by = "ID")

# Annotate neighbors + count
nb_list <- list()

for (sample_name in names(samples)) {
  message("Processing: ", sample_name)
  obj <- samples[[sample_name]]
  
  # compute kNN neighborhood
  obj <- tryCatch(
    compute_neighbors_knn(obj, annotation_col = label_col, k = k_neighbors),
    error = function(e) {
      message("Skipping ", sample_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(obj)) next
  
  # count neighbors
  counts <- tryCatch(
    count_neighbors(obj, label_col = label_col),
    error = function(e) {
      message("Neighbor counting failed for ", sample_name, ": ", e$message)
      NULL
    }
  )
  if (is.null(counts)) next
  
  counts$sample <- sample_name
  nb_list[[sample_name]] <- counts
  
  # store updated sample back 
  samples[[sample_name]] <- obj
}

combined_df <- bind_rows(nb_list)
readr::write_csv(combined_df, out_counts_csv)
message("Wrote: ", out_counts_csv)

# -----------------------------
# Sample summary table (clean version)
# -----------------------------
meta <- RM_ST@meta.data

summary_table <- meta %>%
  group_by(ID) %>%
  summarise(
    `Total # of Spots` = n(),
    `# of Malignant Spots` = sum(grepl("^Cancer_C", .data[[label_col]]), na.rm = TRUE),
    `# of Non-Malignant Spots` = sum(!grepl("^Cancer_C", .data[[label_col]]), na.rm = TRUE),
    `Response` = dplyr::first(Response_Status),
    `Best Response` = dplyr::first(best_overall_response),
    .groups = "drop"
  )

write.csv(summary_table, out_summary_csv, row.names = FALSE)
message("Wrote: ", out_summary_csv)
