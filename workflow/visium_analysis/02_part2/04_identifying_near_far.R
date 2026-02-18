#!/usr/bin/env Rscript
# Part 2.3 — Identify "near" vs "far" spots for the same (N, M) neighbor pairs
#
# Input: list of Seurat objects that already contains:
#   - LABEL_COL (e.g., "Spatial_CT_refined" or "Spatial_CT")
#   - "neighbor_bcs" (comma-separated neighbor barcodes per spot) OR a neighbor graph derived in Part 2.2
#   - one or more *_neighbouring_* columns produced in Part 2.2 (binary or numeric)
#
# Output: same list, with corresponding *_far_* columns added.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(igraph)
  library(tibble)
})

# ============================ USER INPUTS ============================
in_rds  <- "part2/annotated_objects_list_forDEGs.rds"
out_rds <- "part2/annotated_objects_list_forDEGs_withFar.rds"

# Which metadata column defines cell-type / state labels?
LABEL_COL <- "Spatial_CT"  # or "Spatial_CT_refined"

# Graph hop thresholds
HOP_MIN_PRIMARY        <- 3L
HOP_MIN_FALLBACK       <- 2L
MIN_FAR_TARGET         <- 5L
USE_QUANTILE_FALLBACK  <- TRUE
Q_FAR                  <- 0.67   # top one-third most distant among eligible non-neighbors
# ====================================================================

# -------------------------- helpers ---------------------------------
sanitize_label <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}
mk_neigh_col <- function(N, M) paste0(sanitize_label(N), "_neighbouring_", sanitize_label(M))
mk_far_col   <- function(N, M) paste0(sanitize_label(N), "_far_",          sanitize_label(M))

split_clean <- function(x) {
  xs <- strsplit(as.character(x), ",", fixed = TRUE)
  lapply(xs, function(v) v[nzchar(trimws(v))])
}

build_graph <- function(md, neighbor_col = "neighbor_bcs") {
  stopifnot(neighbor_col %in% colnames(md))
  rn <- rownames(md)
  nb_lists <- split_clean(md[[neighbor_col]])
  
  edges <- Map(
    function(src, nbrs) if (length(nbrs)) cbind(src, intersect(nbrs, rn)) else NULL,
    rn, nb_lists
  ) |> do.call(what = rbind)
  
  if (is.null(edges) || nrow(edges) == 0) {
    return(make_empty_graph(n = length(rn), directed = FALSE, names = rn))
  }
  
  graph_from_data_frame(
    data.frame(from = edges[, 1], to = edges[, 2], stringsAsFactors = FALSE),
    directed = FALSE,
    vertices = rn
  )
}

hop_distance_to_class <- function(md, g, label_col, target_label) {
  rn <- rownames(md)
  targets <- rn[md[[label_col]] == target_label]
  out <- rep(Inf, length(rn)); names(out) <- rn
  if (!length(targets)) return(out)
  
  d <- distances(g, v = V(g), to = V(g)[targets], mode = "all")
  out <- if (is.null(dim(d))) as.numeric(d) else apply(d, 1, min, na.rm = TRUE)
  out[is.infinite(out)] <- Inf
  names(out) <- V(g)$name
  out
}

as_logical01 <- function(x) {
  # robust conversion: supports 0/1 numeric, logical, and character "0"/"1"
  if (is.logical(x)) return(x)
  suppressWarnings({
    xn <- as.numeric(as.character(x))
    if (!all(is.na(xn))) return(xn > 0)
  })
  as.logical(x)
}

# -------------------------- main ------------------------------------
objs_in <- readRDS(in_rds)
stopifnot(is.list(objs_in), length(objs_in) > 0)

objs_out <- purrr::imap(objs_in, function(obj, nm) {
  md <- obj@meta.data
  
  if (!(LABEL_COL %in% colnames(md))) {
    message("Skipping ", nm, " (missing LABEL_COL: ", LABEL_COL, ")")
    return(obj)
  }
  if (!("neighbor_bcs" %in% colnames(md))) {
    message("Skipping ", nm, " (missing neighbor_bcs)")
    return(obj)
  }
  
  neigh_cols <- grep("_neighbouring_", colnames(md), value = TRUE, fixed = TRUE)
  if (!length(neigh_cols)) return(obj)
  
  g <- build_graph(md, neighbor_col = "neighbor_bcs")
  
  # parse sanitized pair names from columns
  pairs <- tidyr::separate_wider_delim(
    tibble::tibble(col = neigh_cols),
    col, delim = "_neighbouring_", names = c("N_san", "M_san")
  ) %>%
    mutate(N_san = trimws(N_san), M_san = trimws(M_san))
  
  # map sanitized -> original label, using exact sanitize() matching
  labels_here <- sort(unique(trimws(as.character(md[[LABEL_COL]]))))
  san_map <- setNames(labels_here, sanitize_label(labels_here))
  
  for (i in seq_len(nrow(pairs))) {
    N_san <- pairs$N_san[i]
    M_san <- pairs$M_san[i]
    
    if (!(N_san %in% names(san_map) && M_san %in% names(san_map))) next
    
    N <- san_map[[N_san]]
    M <- san_map[[M_san]]
    
    col_nei <- mk_neigh_col(N, M)
    col_far <- mk_far_col(N, M)
    
    if (!(col_nei %in% colnames(md))) next
    if (col_far %in% colnames(md)) next  # already exists
    
    is_N <- md[[LABEL_COL]] == N
    is_nei <- as_logical01(md[[col_nei]])
    
    # eligible for FAR: N spots that are NOT near M
    eligible <- which(is_N & !is_nei)
    far_flag <- rep(FALSE, nrow(md))
    
    # compute hop distance to malignant class M
    hop <- hop_distance_to_class(md, g, LABEL_COL, M)
    
    # primary threshold
    far_idx <- eligible[hop[eligible] >= HOP_MIN_PRIMARY]
    if (length(far_idx) < MIN_FAR_TARGET) {
      # fallback threshold
      far_idx <- eligible[hop[eligible] >= HOP_MIN_FALLBACK]
    }
    
    if (length(far_idx) < MIN_FAR_TARGET && USE_QUANTILE_FALLBACK) {
      hop_elig <- hop[eligible]
      hop_elig <- hop_elig[is.finite(hop_elig)]
      if (length(hop_elig) > 0) {
        thr <- as.numeric(stats::quantile(hop_elig, probs = Q_FAR, names = FALSE, na.rm = TRUE))
        far_idx <- eligible[hop[eligible] >= thr]
      }
    }
    
    far_flag[far_idx] <- TRUE
    md[[col_far]] <- as.integer(far_flag)  # store as 0/1 for robustness
  }
  
  obj@meta.data <- md
  obj
})

saveRDS(objs_out, out_rds)
message("✅ Wrote: ", out_rds)
