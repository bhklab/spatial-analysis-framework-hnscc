
# ============================================================
# FINAL (NO LEGENDS): ONE CIRCOS PER GROUP (Publication ready)
# - Robust to CellChat versions (ligand/receptor vs interaction_name)
# - Uses LR-level edges (slot.name="net") for top LR per pathway
# - Better label handling (AUTO vs ALL), tighter spacing
# - Smaller arrows, cleaner aesthetics
# - Optional top-K LR labels on chords (OFF by default)
# ============================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(circlize)
  library(grid)
})

# ----------------------------
# USER SETTINGS
# ----------------------------
TOP_PER_PATHWAY <- 5
DPI   <- 900
FIG_W <- 6.6
FIG_H <- 6.6

out_dir <- file.path("figures","fig4_cellchat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Inputs (edit these paths)
# ----------------------------
cellchat_tier2_rds  <- "/part3/cellchat_tier2.rds"
cellchat_immune_rds <- "/part3/cellchat_tier2c_immune.rds"

stopifnot(file.exists(cellchat_tier2_rds), file.exists(cellchat_immune_rds))
cellchat_tier2  <- readRDS(cellchat_tier2_rds)
cellchat_immune <- readRDS(cellchat_immune_rds)


safe_name <- function(x) gsub("[^A-Za-z0-9]+", "", x)

# ============================================================
# Helper: robustly standardize CellChat communication table
# - Prefer slot.name="net" (LR-level). netP can lose ligand/receptor columns.
# - If ligand/receptor columns absent, parse interaction_name / interaction_name_2
# ============================================================
.standardize_comm_df <- function(comm) {
  cn <- colnames(comm)
  
  # must-have: source/target
  if (!all(c("source","target") %in% cn)) {
    stop("subsetCommunication output missing 'source'/'target'. Columns are: ",
         paste(cn, collapse = ", "))
  }
  
  # pathway column
  path_col <- intersect(cn, c("pathway_name","pathway"))[1]
  if (is.na(path_col)) {
    stop("Cannot find pathway column in subsetCommunication output. Columns are: ",
         paste(cn, collapse = ", "))
  }
  
  # score column: prob or weight (prob is typical)
  score_col <- intersect(cn, c("prob","weight","prob.mean","prob_mean"))[1]
  if (is.na(score_col)) {
    stop("Cannot find score column (prob/weight/prob.mean) in subsetCommunication output. Columns are: ",
         paste(cn, collapse = ", "))
  }
  
  # ligand/receptor columns (preferred)
  lig_col <- intersect(cn, c("ligand","ligand.gene","ligand_gene","ligand_symbol","ligand_name"))[1]
  rec_col <- intersect(cn, c("receptor","receptor.gene","receptor_gene","receptor_symbol","receptor_name"))[1]
  
  # fallback: interaction_name
  int_col <- intersect(cn, c("interaction_name","interaction_name_2","interaction"))[1]
  
  if (!is.na(lig_col) && !is.na(rec_col)) {
    out <- comm %>%
      transmute(
        source       = .data[["source"]],
        target       = .data[["target"]],
        pathway_name = .data[[path_col]],
        ligand_gene  = .data[[lig_col]],
        receptor_gene= .data[[rec_col]],
        score        = .data[[score_col]],
        pval         = if ("pval" %in% cn) .data[["pval"]] else NA_real_,
        p_adjust     = if ("p.adjust" %in% cn) .data[["p.adjust"]] else NA_real_
      )
    return(out)
  }
  
  if (!is.na(int_col)) {
    # Typical CellChat interaction_name format is "LIGAND_RECEPTOR"
    # If yours uses "-" or other delimiters, this still tries reasonably.
    tmp <- comm[[int_col]]
    tmp <- as.character(tmp)
    
    split_lr <- strsplit(tmp, "_", fixed = TRUE)
    
    lig <- vapply(split_lr, function(z) if (length(z) >= 1) z[[1]] else NA_character_, character(1))
    rec <- vapply(split_lr, function(z) if (length(z) >= 2) paste(z[-1], collapse = "_") else NA_character_, character(1))
    
    # If underscore split fails (no receptor), try "-" split
    bad <- is.na(rec) | rec == ""
    if (any(bad)) {
      split_lr2 <- strsplit(tmp[bad], "-", fixed = TRUE)
      lig[bad] <- vapply(split_lr2, function(z) if (length(z) >= 1) z[[1]] else NA_character_, character(1))
      rec[bad] <- vapply(split_lr2, function(z) if (length(z) >= 2) paste(z[-1], collapse = "-") else NA_character_, character(1))
    }
    
    out <- comm %>%
      transmute(
        source       = .data[["source"]],
        target       = .data[["target"]],
        pathway_name = .data[[path_col]],
        ligand_gene  = lig,
        receptor_gene= rec,
        score        = .data[[score_col]],
        pval         = if ("pval" %in% cn) .data[["pval"]] else NA_real_,
        p_adjust     = if ("p.adjust" %in% cn) .data[["p.adjust"]] else NA_real_
      )
    
    return(out)
  }
  
  stop(
    "Could not find ligand/receptor columns OR interaction_name column.\n",
    "Columns are: ", paste(cn, collapse = ", "),
    "\nTip: run `colnames(subsetCommunication(cellchat_obj, slot.name='net'))` and check what columns exist."
  )
}

# ============================================================
# 1) Extract observed LR edges for selected pathways & direction
#    then keep top N per pathway (by summed prob/weight)
# ============================================================
get_observed_topN_per_pathway <- function(cellchat_obj,
                                          pathways,
                                          sources = NULL,
                                          targets = NULL,
                                          top_n = 5,
                                          pval_cutoff = NULL,
                                          slot.name = "net") {
  
  # IMPORTANT: use net for LR-level edges (best for "top LR per pathway")
  comm_raw <- subsetCommunication(cellchat_obj, slot.name = slot.name)
  
  df0 <- .standardize_comm_df(comm_raw)
  
  df <- df0 %>% filter(pathway_name %in% pathways)
  if (!is.null(sources)) df <- df %>% filter(source %in% sources)
  if (!is.null(targets)) df <- df %>% filter(target %in% targets)
  
  if (!is.null(pval_cutoff)) {
    df <- df %>%
      filter(is.na(pval) | pval <= pval_cutoff) %>%
      filter(is.na(p_adjust) | p_adjust <= pval_cutoff)
  }
  
  df2 <- df %>%
    filter(!is.na(ligand_gene), !is.na(receptor_gene), !is.na(score)) %>%
    mutate(LR = paste0(ligand_gene, "-", receptor_gene))
  
  top_lr <- df2 %>%
    group_by(pathway_name, ligand_gene, receptor_gene, LR) %>%
    summarise(score = sum(score, na.rm = TRUE), .groups = "drop") %>%
    group_by(pathway_name) %>%
    arrange(desc(score), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  summary <- top_lr %>%
    group_by(pathway_name) %>%
    summarise(
      n_LR = n(),
      total_prob_topN = sum(score, na.rm = TRUE),
      LR_list = paste(LR, collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(match(pathway_name, pathways))
  
  list(top_lr_long = top_lr, summary = summary)
}

# ============================================================
# 2) Publication-ready plot (NO LEGENDS)
#    - Better label strategy
#    - Smaller arrows
#    - Cleaner spacing
#    - Optional chord LR labels (top K overall)
# ============================================================
plot_group_circos_topLR_pub_nolegend <- function(top_lr_long,
                                                 pathways_order,
                                                 title_text = NULL,
                                                 ligand_col = "#CC79A7",
                                                 receptor_col = "#D55E00",
                                                 pathway_cols = NULL,
                                                 start_degree = 90,
                                                 small_gap = 1.0,
                                                 big_gap = 10,
                                                 track_height = 0.10,
                                                 label_cex = 0.55,
                                                 title_cex = 20,
                                                 chord_alpha = 0.30,
                                                 chord_lwd = 1.05,
                                                 weight_transform = c("none","sqrt","log1p"),
                                                 # label strategy
                                                 label_mode = c("auto","all"),
                                                 label_top_n = 18,
                                                 label_min_frac = 0.015,
                                                 # arrows (smaller)
                                                 arrow_diffHeight = 0.25,
                                                 arrow_length     = 0.18,
                                                 # label distance from chords (mm)
                                                 label_inset_mm = 0.06,
                                                 # Optional: chord LR labels
                                                 show_chord_labels = FALSE,
                                                 chord_label_top_k = 10,
                                                 chord_label_cex   = 0.45) {
  
  weight_transform <- match.arg(weight_transform)
  label_mode <- match.arg(label_mode)
  
  df <- top_lr_long
  if (nrow(df) == 0) {
    message("⚠️ No LR edges to plot.")
    return(invisible(NULL))
  }
  
  df <- df %>% mutate(pathway_name = factor(pathway_name, levels = pathways_order))
  
  # weights -> chord width
  w <- df$score
  if (weight_transform == "sqrt")  w <- sqrt(w)
  if (weight_transform == "log1p") w <- log1p(w)
  
  ligs <- unique(df$ligand_gene)
  recs <- unique(df$receptor_gene)
  
  lig_ids <- paste0("L_", ligs)
  rec_ids <- paste0("R_", recs)
  
  links_df <- df %>%
    mutate(score_plot = w) %>%
    transmute(
      from  = paste0("L_", ligand_gene),
      to    = paste0("R_", receptor_gene),
      value = score_plot,
      pathway_name = as.character(pathway_name),
      LR = paste0(ligand_gene, "-", receptor_gene)
    )
  
  ord <- c(lig_ids, rev(paste0("R_", recs)))
  
  # sector colors
  grid_cols <- c(
    setNames(rep(ligand_col, length(lig_ids)), lig_ids),
    setNames(rep(receptor_col, length(rec_ids)), paste0("R_", recs))
  )
  
  # pathway colors
  if (is.null(pathway_cols)) {
    base_pal <- c("#4E79A7","#F28E2B","#59A14F","#E15759","#B07AA1",
                  "#76B7B2","#9C755F","#EDC948","#AFB0B2","#1F77B4")
    pw <- pathways_order
    pathway_cols <- if (length(pw) <= length(base_pal)) {
      setNames(base_pal[seq_along(pw)], pw)
    } else {
      setNames(colorRampPalette(base_pal)(length(pw)), pw)
    }
  } else {
    stopifnot(all(pathways_order %in% names(pathway_cols)))
    pathway_cols <- pathway_cols[pathways_order]
  }
  
  link_cols <- pathway_cols[links_df$pathway_name]
  link_cols <- grDevices::adjustcolor(link_cols, alpha.f = chord_alpha)
  
  # label selection
  if (label_mode == "all") {
    label_keep <- ord
  } else {
    flow_from <- links_df %>% group_by(from) %>% summarise(v=sum(value), .groups="drop")
    flow_to   <- links_df %>% group_by(to)   %>% summarise(v=sum(value), .groups="drop")
    flow <- bind_rows(flow_from %>% rename(sector=from),
                      flow_to   %>% rename(sector=to)) %>%
      group_by(sector) %>% summarise(v=sum(v), .groups="drop")
    flow$frac <- flow$v / sum(flow$v)
    
    keep_by_frac <- flow$sector[flow$frac >= label_min_frac]
    keep_by_rank <- flow %>% arrange(desc(v)) %>% slice_head(n = min(label_top_n, nrow(flow))) %>% pull(sector)
    label_keep <- unique(c(keep_by_frac, keep_by_rank))
  }
  
  # gaps
  n_l <- length(lig_ids)
  n_r <- length(recs)
  gap_after <- c(
    rep(small_gap, max(0, n_l - 1)), big_gap,
    rep(small_gap, max(0, n_r - 1)), big_gap
  )
  
  # draw
  circos.clear()
  circos.par(
    start.degree = start_degree,
    gap.after    = gap_after,
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim  = c(-1.04, 1.04),
    canvas.ylim  = c(-1.06, 1.06),
    points.overflow.warning = FALSE
  )
  
  chordDiagramFromDataFrame(
    df = links_df %>% select(from, to, value, pathway_name),
    order = ord,
    grid.col = grid_cols,
    transparency = 0,
    col = link_cols,
    directional = 1,
    direction.type = "arrows",
    diffHeight = arrow_diffHeight,
    link.arr.length = arrow_length,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = track_height),
    link.border = NA,
    link.lwd = chord_lwd,
    link.sort = TRUE,
    link.decreasing = TRUE
  )
  
  # sector labels (close to chord anchors)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    nm   <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    if (!nm %in% label_keep) return(NULL)
    
    nm_clean <- sub("^L_", "", nm)
    nm_clean <- sub("^R_", "", nm_clean)
    
    circos.text(
      x = mean(xlim),
      y = ylim[1] + mm_y(label_inset_mm),
      labels = nm_clean,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = label_cex
    )
  }, bg.border = NA)
  
  # Optional: label top-K LR pairs (center-ish). Use sparingly for publication.
  if (isTRUE(show_chord_labels)) {
    lab_df <- links_df %>%
      arrange(desc(value)) %>%
      distinct(from, to, .keep_all = TRUE) %>%
      slice_head(n = min(chord_label_top_k, nrow(.)))
    
    # Place labels slightly inside the circle near the center (global positions).
    # (This is a heuristic; looks good for top few labels.)
    for (i in seq_len(nrow(lab_df))) {
      # approximate mid-angle between sectors
      f <- lab_df$from[i]; t <- lab_df$to[i]
      af <- get.cell.meta.data("cell.start.degree", sector.index = f)
      bf <- get.cell.meta.data("cell.end.degree",   sector.index = f)
      at <- get.cell.meta.data("cell.start.degree", sector.index = t)
      bt <- get.cell.meta.data("cell.end.degree",   sector.index = t)
      ang <- mean(c(mean(c(af,bf)), mean(c(at,bt))))
      rad <- 0.58
      circos.text(
        x = rad * cospi(ang/180),
        y = rad * sinpi(ang/180),
        labels = lab_df$LR[i],
        cex = chord_label_cex,
        facing = "inside",
        niceFacing = TRUE
      )
    }
  }
  
  if (!is.null(title_text)) {
    grid.text(title_text, x = unit(0.5, "npc"), y = unit(0.97, "npc"),
              gp = gpar(fontsize = title_cex, fontface = "bold"))
  }
  
  invisible(NULL)
}

# ============================================================
# 3) Wrapper: one group -> one CSV + one PNG + one PDF
# ============================================================
run_one_group <- function(group_label,
                          cellchat_obj,
                          pathways,
                          sources,
                          targets,
                          ligand_col,
                          receptor_col,
                          top_per_pathway = 5,
                          weight_transform = "sqrt",
                          label_inset_mm = 0.06,
                          label_mode = "auto",
                          label_top_n = 18,
                          label_min_frac = 0.015,
                          # arrows
                          arrow_diffHeight = 0.25,
                          arrow_length     = 0.18,
                          # optional chord labels
                          show_chord_labels = FALSE,
                          chord_label_top_k = 10) {
  
  message("\n=== ", group_label, " ===")
  
  res <- get_observed_topN_per_pathway(
    cellchat_obj = cellchat_obj,
    pathways     = pathways,
    sources      = sources,
    targets      = targets,
    top_n        = top_per_pathway,
    slot.name    = "net"   # << keep LR-level table
  )
  
  # CSV summary
  csv_path <- file.path(out_dir, paste0("Observed_", safe_name(group_label),
                                        "_Top", top_per_pathway, "PerPathway.csv"))
  write.csv(res$summary, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  message("✅ Saved CSV: ", csv_path)
  
  # PNG
  png(file.path(out_dir, paste0("Circos_", safe_name(group_label),
                                "_Top", top_per_pathway, "PerPathway.png")),
      width = FIG_W, height = FIG_H, units = "in", res = DPI, bg = "white")
  par(mar = c(0, 0, 0, 0))
  plot_group_circos_topLR_pub_nolegend(
    top_lr_long      = res$top_lr_long,
    pathways_order   = pathways,
    title_text       = paste0(group_label, " | Top ", top_per_pathway, " LR / pathway"),
    ligand_col       = ligand_col,
    receptor_col     = receptor_col,
    weight_transform = weight_transform,
    label_mode       = label_mode,
    label_top_n      = label_top_n,
    label_min_frac   = label_min_frac,
    label_cex        = 0.55,
    label_inset_mm   = label_inset_mm,
    chord_alpha      = 0.30,
    chord_lwd        = 1.05,
    title_cex        = 20,
    track_height     = 0.10,
    arrow_diffHeight = arrow_diffHeight,
    arrow_length     = arrow_length,
    show_chord_labels= show_chord_labels,
    chord_label_top_k= chord_label_top_k
  )
  dev.off()
  
  # PDF
  pdf(file.path(out_dir, paste0("Circos_", safe_name(group_label),
                                "_Top", top_per_pathway, "PerPathway.pdf")),
      width = FIG_W, height = FIG_H)
  par(mar = c(0, 0, 0, 0))
  plot_group_circos_topLR_pub_nolegend(
    top_lr_long      = res$top_lr_long,
    pathways_order   = pathways,
    title_text       = NULL,
    ligand_col       = ligand_col,
    receptor_col     = receptor_col,
    weight_transform = weight_transform,
    label_mode       = label_mode,
    label_top_n      = label_top_n,
    label_min_frac   = label_min_frac,
    label_cex        = 0.55,
    label_inset_mm   = label_inset_mm,
    chord_alpha      = 0.30,
    chord_lwd        = 1.05,
    title_cex        = 20,
    track_height     = 0.10,
    arrow_diffHeight = arrow_diffHeight,
    arrow_length     = arrow_length,
    show_chord_labels= show_chord_labels,
    chord_label_top_k= chord_label_top_k
  )
  dev.off()
  
  message("✅ Saved PNG/PDF: ", group_label)
  invisible(res)
}

# ============================================================
# DEFINE GROUPS (must match your idents)
# ============================================================
caf_groups   <- c("non-activated fibroblast","ecm_myCAF","IFNg-iCAF","IL-iCAF","acto-myCAF")
malig_groups <- c("cycling","TC","neutrophil_inflamed","fibrovascular","pEMT")

# ============================================================
# RUN GROUPS (Publication-ready defaults)
# - label_mode="auto" is clean; use "all" if you want every gene labeled
# - arrows are smaller by default now
# ============================================================

res_cafs_to_cycling <- run_one_group(
  group_label  = "CAFs->cycling",
  cellchat_obj = cellchat_tier2,
  pathways     = c("NGF","PERIOSTIN","COLLAGEN","LAMININ","FN1"),
  sources      = caf_groups,
  targets      = "cycling",
  ligand_col   = "#CC79A7",
  receptor_col = "#D55E00",
  top_per_pathway  = TOP_PER_PATHWAY,
  weight_transform = "sqrt",
  label_mode       = "auto",   # change to "all" if you want all gene labels
  label_inset_mm   = 0.06,
  arrow_length     = 0.18,
  arrow_diffHeight = 0.25,
  show_chord_labels = FALSE
)

res_cycling_to_cafs <- run_one_group(
  group_label  = "cycling->CAFs",
  cellchat_obj = cellchat_tier2,
  pathways     = c("KIT","AGRN","WNT","IFN-I","LAMININ"),
  sources      = "cycling",
  targets      = caf_groups,
  ligand_col   = "#D55E00",
  receptor_col = "#CC79A7",
  top_per_pathway  = TOP_PER_PATHWAY,
  weight_transform = "sqrt",
  label_mode       = "auto",
  label_inset_mm   = 0.06,
  arrow_length     = 0.18,
  arrow_diffHeight = 0.25,
  show_chord_labels = FALSE
)

res_ifntam_to_malig <- run_one_group(
  group_label  = "IFN-TAM->malignant",
  cellchat_obj = cellchat_immune,
  pathways     = c("BAFF","FN1","GALECTIN"),
  sources      = "IFN-TAM",
  targets      = malig_groups,
  ligand_col   = "#8B4C9A",
  receptor_col = "#D55E00",
  top_per_pathway  = TOP_PER_PATHWAY,
  weight_transform = "sqrt",
  label_mode       = "auto",
  label_inset_mm   = 0.06,
  arrow_length     = 0.18,
  arrow_diffHeight = 0.25,
  show_chord_labels = FALSE
)

# IFN-TAM -> pooled TAM subtypes (excluding IFN-TAM)
all_idents_immune <- levels(cellchat_immune@idents)
tam_subtypes_other <- all_idents_immune[grepl("TAM", all_idents_immune, ignore.case = TRUE)]
tam_subtypes_other <- setdiff(tam_subtypes_other, c("IFN-TAM","IFN.TAM","IFN_TAM","IFNTAM"))

if (length(tam_subtypes_other) == 0) {
  warning("No TAM targets detected for pooled TAM subtypes. Check cellchat_immune@idents labels.")
} else {
  message("Pooled TAM targets (excluding IFN-TAM): ", paste(tam_subtypes_other, collapse = ", "))
  
  res_ifntam_to_tams <- run_one_group(
    group_label  = "IFN-TAM->TAMsOther",
    cellchat_obj = cellchat_immune,
    pathways     = c("BAFF","TNF","FN1","TIGIT","GALECTIN"),
    sources      = "IFN-TAM",
    targets      = tam_subtypes_other,
    ligand_col   = "#8B4C9A",
    receptor_col = "#7A7A7A",
    top_per_pathway  = TOP_PER_PATHWAY,
    weight_transform = "sqrt",
    label_mode       = "auto",
    label_inset_mm   = 0.06,
    arrow_length     = 0.18,
    arrow_diffHeight = 0.25,
    show_chord_labels = FALSE
  )
}

