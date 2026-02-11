#!/usr/bin/env Rscript
# Part 2.4 — Differential expression: Near vs Far (within N) conditioned on malignant state M
#
# Requires Part 2.3 output list that includes *_neighbouring_* and *_far_* columns.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(readr)
})

# ============================ USER INPUTS ============================
in_rds  <- "part2/annotated_objects_list_forDEGs_withFar.rds"
out_dir <- "part2/2.4_DEG_near_vs_far"

# Which label column defines cell types/states?
LABEL_COL <- "Spatial_CT"  # or "Spatial_CT_refined"

# Pair list (N = neighbor cell type, M = malignant state)
combos <- tribble(
  ~N, ~M
  # Example:
  # "mCAFs", "Cancer_C1"
)

# Optional filter to run only PRE objects (by list name suffix or metadata)
FILTER_MODE <- c("none", "name_suffix", "metadata")[1]
PRE_SUFFIX  <- "_Pre$"          # used if FILTER_MODE == "name_suffix"
META_FIELD  <- "Pre_Post"       # used if FILTER_MODE == "metadata"
META_VALUE  <- "Pre"

# DE parameters
ASSAY_NAME       <- "Spatial"   # your data layer for DE
SLOT_NAME        <- "data"      # typically "data" (log-normalized) for FindMarkers
MIN_PCT          <- 0.25
LOGFC_MIN        <- 0.25
PADJ_MAX         <- 0.05
LOGFC_MIN_REPORT <- 0.58        # for reporting/overlap summaries (~1.5x)
MIN_CELLS_GROUP  <- 10
TEST_USE         <- "wilcox"    # publication-safe default
# ====================================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- helpers ---------------------------------
sanitize_label <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("_+", "_", y)
  y <- sub("^_", "", y); y <- sub("_$", "", y)
  y
}
mk_neigh_col <- function(N, M) paste0(sanitize_label(N), "_neighbouring_", sanitize_label(M))
mk_far_col   <- function(N, M) paste0(sanitize_label(N), "_far_",          sanitize_label(M))
as_logical01 <- function(x) {
  if (is.logical(x)) return(x)
  suppressWarnings({
    xn <- as.numeric(as.character(x))
    if (!all(is.na(xn))) return(xn > 0)
  })
  as.logical(x)
}

filter_objects <- function(objs) {
  if (FILTER_MODE == "none") return(objs)
  
  if (FILTER_MODE == "name_suffix") {
    keep <- grepl(PRE_SUFFIX, names(objs))
    return(objs[keep])
  }
  
  if (FILTER_MODE == "metadata") {
    keep <- purrr::imap_lgl(objs, function(obj, nm) {
      md <- obj@meta.data
      if (!(META_FIELD %in% colnames(md))) return(FALSE)
      as.character(md[[META_FIELD]][1]) == META_VALUE
    })
    return(objs[keep])
  }
  
  objs
}

# -------------------------- load ------------------------------------
objs_all <- readRDS(in_rds)
stopifnot(is.list(objs_all), length(objs_all) > 0)

objs <- filter_objects(objs_all)
stopifnot(length(objs) > 0)

# -------------------------- per-sample DE ---------------------------
per_sample_pair_results <- list()

for (k in seq_len(nrow(combos))) {
  N <- combos$N[k]
  M <- combos$M[k]
  
  col_nei <- mk_neigh_col(N, M)
  col_far <- mk_far_col(N, M)
  
  pair_key <- paste0(sanitize_label(N), "__Neighbor_vs_Far__", sanitize_label(M))
  message("\n=== Pair: ", N, " (Neighbor vs Far), conditioned on M=", M, " ===")
  
  pair_dir <- file.path(out_dir, pair_key)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  
  pair_res <- purrr::imap(objs, function(obj, nm) {
    md <- obj@meta.data
    
    # checks
    if (!(LABEL_COL %in% colnames(md))) return(NULL)
    if (!all(c(col_nei, col_far) %in% colnames(md))) return(NULL)
    
    is_N   <- md[[LABEL_COL]] == N
    is_nei <- as_logical01(md[[col_nei]])
    is_far <- as_logical01(md[[col_far]])
    
    cells_nei <- rownames(md)[is_N & is_nei]
    cells_far <- rownames(md)[is_N & is_far]
    
    if (length(cells_nei) < MIN_CELLS_GROUP || length(cells_far) < MIN_CELLS_GROUP) {
      return(NULL)
    }
    
    DefaultAssay(obj) <- ASSAY_NAME
    
    # Run DE: Neighbor (ident.1) vs Far (ident.2)
    deg <- FindMarkers(
      object = obj,
      ident.1 = cells_nei,
      ident.2 = cells_far,
      assay   = ASSAY_NAME,
      slot    = SLOT_NAME,
      min.pct = MIN_PCT,
      logfc.threshold = LOGFC_MIN,
      test.use = TEST_USE
    ) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(
        sample = nm,
        N = N,
        M = M,
        n_neighbor = length(cells_nei),
        n_far      = length(cells_far)
      )
    
    # write per-sample DE
    out_csv <- file.path(pair_dir, paste0("DEG_", sanitize_label(nm), ".csv"))
    readr::write_csv(deg, out_csv)
    
    deg
  })
  
  pair_res <- purrr::compact(pair_res)
  per_sample_pair_results[[pair_key]] <- pair_res
}

# -------------------------- overlap summaries ------------------------
overlap_summaries <- list()

for (pair_key in names(per_sample_pair_results)) {
  res_list <- per_sample_pair_results[[pair_key]]
  if (length(res_list) == 0) next
  
  # bind all sample results for this pair
  all_deg <- bind_rows(res_list)
  
  # choose column name depending on Seurat version
  lfc_col <- if ("avg_log2FC" %in% colnames(all_deg)) "avg_log2FC" else "avg_logFC"
  
  sig <- all_deg %>%
    filter(!is.na(p_val_adj), p_val_adj <= PADJ_MAX, abs(.data[[lfc_col]]) >= LOGFC_MIN_REPORT)
  
  if (nrow(sig) == 0) next
  
  # count recurrence per gene, direction
  occ <- sig %>%
    mutate(direction = ifelse(.data[[lfc_col]] > 0, "Up_in_Neighbor", "Up_in_Far")) %>%
    count(gene, direction) %>%
    tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
    mutate(
      total_sig = Up_in_Neighbor + Up_in_Far,
      consensus = case_when(
        Up_in_Neighbor > 0 & Up_in_Far == 0 ~ "Neighbor_only",
        Up_in_Far > 0 & Up_in_Neighbor == 0 ~ "Far_only",
        Up_in_Neighbor > 0 & Up_in_Far > 0 ~ "Mixed",
        TRUE ~ "NA"
      )
    ) %>%
    arrange(desc(total_sig), desc(Up_in_Neighbor), desc(Up_in_Far))
  
  out_csv <- file.path(out_dir, paste0("OVERLAP_", pair_key, ".csv"))
  readr::write_csv(occ, out_csv)
  overlap_summaries[[pair_key]] <- occ
  message("Wrote overlap summary: ", out_csv)
}

# Combine all overlap summaries
if (length(overlap_summaries) > 0) {
  all_pairs_summary <- purrr::imap_dfr(overlap_summaries, function(df, key) {
    df %>% mutate(pair = key) %>%
      select(pair, gene, Up_in_Neighbor, Up_in_Far, total_sig, consensus)
  })
  readr::write_csv(all_pairs_summary, file.path(out_dir, "OVERLAP_all_pairs_summary.csv"))
}

# Save RDS bundle
saveRDS(
  list(
    inputs = list(in_rds = in_rds, combos = combos, LABEL_COL = LABEL_COL),
    per_sample_pair_results = per_sample_pair_results,
    overlap_summaries = overlap_summaries,
    params = list(
      ASSAY_NAME = ASSAY_NAME, SLOT_NAME = SLOT_NAME,
      MIN_PCT = MIN_PCT, LOGFC_MIN = LOGFC_MIN, PADJ_MAX = PADJ_MAX,
      LOGFC_MIN_REPORT = LOGFC_MIN_REPORT, MIN_CELLS_GROUP = MIN_CELLS_GROUP,
      FILTER_MODE = FILTER_MODE, PRE_SUFFIX = PRE_SUFFIX,
      META_FIELD = META_FIELD, META_VALUE = META_VALUE,
      TEST_USE = TEST_USE
    )
  ),
  file.path(out_dir, "DEG_near_vs_far_results.rds")
)

message("✅ Done. Results in: ", out_dir)
