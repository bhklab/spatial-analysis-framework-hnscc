# ============================================================
# Supplementary Figure 1
# Publication polished script (paths parameterized, no setwd)
# ============================================================


# ---- shared helpers ----
# Assumes this script lives in the same folder as 00_common.R.
# If not, change the path below (e.g., source("code/04_part4/00_common.R")).
source("00_common.R")
theme_set(theme_pub())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
  library(tidyr)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------
input_rds_1 <- "/part2/RM_ST.rds"

# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "suppl_fig1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


#Pie Chart
###############################################
## Pie charts per Sample_ID_renamed (Baseline-only vs Paired)
## - 4 broad categories: Endothelial / Fibroblast / Immune / Malignant
## - Auto-detect baseline label (BL vs Pre)
## - Faceted pies with percent labels
## - Publication style (Cancer Discovery)
###############################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(patchwork)
  library(tidyr)
})
# ------------------------------------------------------------
# 0) RM_ST assumed already loaded
# RM_ST <- readRDS(input_rds_1)
# ------------------------------------------------------------

# ----------------------------
# 0.1) Ensure Sample_ID_renamed exists (skip if already in your object)
# ----------------------------
if (!"Sample_ID_renamed" %in% colnames(RM_ST@meta.data)) {
  message("Sample_ID_renamed not found. Creating it using Post->OnT logic...")
  
  keep_post <- c("LIB15_Post", "LIB17_Post")
  old_ids <- RM_ST$ID
  new_ids <- old_ids
  new_ids <- gsub("_Post$", "_OnT", new_ids)
  new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]
  RM_ST$Sample_ID_renamed <- new_ids
}

# ----------------------------
# 1) Create 4 broad categories
# ----------------------------
if (!"Spatial_CT_broad" %in% colnames(RM_ST@meta.data)) {
  stop("Spatial_CT_broad not found. Please create it first (as in your earlier script).")
}

RM_ST$CellType4 <- dplyr::case_when(
  RM_ST$Spatial_CT_broad %in% c("Malignant cells")    ~ "Malignant",
  RM_ST$Spatial_CT_broad %in% c("Fibroblasts")        ~ "Fibroblast",
  RM_ST$Spatial_CT_broad %in% c("Endothelial cells")  ~ "Endothelial",
  TRUE                                                ~ "Immune"
)

RM_ST$CellType4 <- factor(
  RM_ST$CellType4,
  levels = c("Endothelial", "Fibroblast", "Immune", "Malignant")
)

# ----------------------------
# 2) Parse patient + timepoint from Sample_ID_renamed
# ----------------------------
meta_df <- RM_ST@meta.data %>%
  mutate(
    Sample_ID_renamed = as.character(Sample_ID_renamed),
    
    # Timepoint = suffix after last underscore (BL, Pre, OnT, Post, C2T, C3T, etc.)
    Timepoint = str_extract(Sample_ID_renamed, "(?<=_)[^_]+$"),
    
    # Patient/root = everything before last underscore
    Patient   = str_replace(Sample_ID_renamed, "_[^_]+$", ""),
    
    # Two-line label like example
    Sample_label = paste0(Patient, "\n", Timepoint)
  )

# ----------------------------
# 2.5) Auto-detect baseline label ("BL" vs "Pre")
# ----------------------------
baseline_tp <- dplyr::case_when(
  "BL"  %in% meta_df$Timepoint ~ "BL",
  "Pre" %in% meta_df$Timepoint ~ "Pre",
  TRUE ~ NA_character_
)

if (is.na(baseline_tp)) {
  stop(
    "Could not find baseline timepoint label. Expected 'BL' or 'Pre'. Found: ",
    paste(sort(unique(meta_df$Timepoint)), collapse = ", ")
  )
}
message("Using baseline timepoint = ", baseline_tp)


# 2.6) Define on-treatment and post-treatment labels
# ----------------------------
on_candidates <- c("OnT", "C2T", "C3T")
on_tps <- intersect(on_candidates, unique(meta_df$Timepoint))

post_tp <- if ("Post" %in% unique(meta_df$Timepoint)) "Post" else NA_character_

message("Detected on-treatment timepoints: ", ifelse(length(on_tps)==0, "NONE", paste(on_tps, collapse=", ")))
message("Detected post-treatment timepoint: ", ifelse(is.na(post_tp), "NONE", post_tp))



# ----------------------------
# 3) Patient grouping into 3 columns (UPDATED PAIRING)
# ----------------------------
patient_tp <- meta_df %>%
  distinct(Patient, Timepoint) %>%
  group_by(Patient) %>%
  summarise(
    has_pre  = any(Timepoint == baseline_tp),
    has_on   = any(Timepoint %in% on_tps),
    has_post = ifelse(is.na(post_tp), FALSE, any(Timepoint == post_tp)),
    n_tp     = n(),
    .groups = "drop"
  )

# Column 1: Pre-treatment only (baseline only)
pre_only_patients <- patient_tp %>%
  filter(has_pre, !has_on, !has_post, n_tp == 1) %>%
  pull(Patient)

# Column 2: Paired (Pre + On) OR (Pre + Post) OR (Pre + On + Post)
paired_patients <- patient_tp %>%
  filter(has_pre, (has_on | has_post)) %>%
  pull(Patient)

# Column 3: Post-treatment only (no baseline, only Post)
post_only_patients <- patient_tp %>%
  filter(!has_pre, has_post, n_tp == 1) %>%
  pull(Patient)

# ----------------------------
# 4) Build per-sample proportions (one pie per Sample_ID_renamed)
# ----------------------------
pie_df <- meta_df %>%
  group_by(Patient, Timepoint, Sample_ID_renamed, Sample_label, CellType4) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample_ID_renamed) %>%
  mutate(
    pct = n / sum(n),
    pct_lab = ifelse(pct >= 0.05, percent(pct, accuracy = 1), "")  # label >= 5%
  ) %>%
  ungroup()

# ----------------------------
# 5) Ordering
# ----------------------------
tp_levels <- unique(c(baseline_tp, on_candidates, "Post"))
all_tp <- unique(pie_df$Timepoint)
tp_levels2 <- unique(c(tp_levels, setdiff(sort(all_tp), tp_levels)))

pie_df <- pie_df %>%
  mutate(
    Timepoint = factor(Timepoint, levels = tp_levels2),
    Patient   = factor(Patient, levels = sort(unique(Patient)))
  )

pre_only_order <- pie_df %>%
  filter(Patient %in% pre_only_patients, Timepoint == baseline_tp) %>%
  distinct(Sample_label, Patient, Timepoint) %>%
  arrange(Patient) %>%
  pull(Sample_label)

# Paired: show Pre + (On if present) + (Post if present)
paired_order <- pie_df %>%
  filter(Patient %in% paired_patients,
         Timepoint %in% c(baseline_tp, on_tps, post_tp)) %>%
  distinct(Sample_label, Patient, Timepoint) %>%
  arrange(Patient, Timepoint) %>%
  pull(Sample_label)

post_only_order <- pie_df %>%
  filter(Patient %in% post_only_patients, Timepoint == post_tp) %>%
  distinct(Sample_label, Patient, Timepoint) %>%
  arrange(Patient) %>%
  pull(Sample_label)

# ----------------------------
# 6) Colors (journal-friendly)
# ----------------------------
celltype4_cols <- c(
  "Endothelial" = "#D9EF8B",  # light yellow-green
  "Fibroblast"  = "#C8BE9A",  # warm khaki/tan
  "Immune"      = "#2C6DB2",  # deep blue
  "Malignant"   = "#E67C73"   # salmon/red
)


# ----------------------------
# 7) Fixed facet geometry so pies are SAME SIZE across panels
# ----------------------------
NROW_FIXED <- 4

n_fac_pre  <- length(pre_only_order)
n_fac_pair <- length(paired_order)
n_fac_post <- length(post_only_order)

NCOL_pre  <- max(1, ceiling(n_fac_pre  / NROW_FIXED))
NCOL_pair <- max(1, ceiling(n_fac_pair / NROW_FIXED))
NCOL_post <- max(1, ceiling(n_fac_post / NROW_FIXED))

NCOL_global <- max(NCOL_pre, NCOL_pair, NCOL_post)

message("Facet grid (nrow fixed = ", NROW_FIXED, "): ",
        "Pre-only ncol=", NCOL_pre, "; Paired ncol=", NCOL_pair, "; Post-only ncol=", NCOL_post)


# ----------------------------
# 8) Helper: pad with dummy facets to fill a grid (critical for Post-only)
# ----------------------------
pad_to_grid <- function(df, label_levels, nrow_fixed, ncol_fixed) {
  target_n <- nrow_fixed * ncol_fixed
  have_n   <- length(label_levels)
  if (have_n >= target_n) return(list(df = df, levels = label_levels))
  
  add_n <- target_n - have_n
  dummy_levels <- paste0(".__dummy__", seq_len(add_n))
  
  dummy_df <- expand_grid(
    Sample_label = dummy_levels,
    CellType4    = levels(df$CellType4)
  ) %>%
    mutate(
      pct = 0,
      pct_lab = "",
      lab_col = "black"
    )
  
  out_df <- bind_rows(df, dummy_df)
  out_levels <- c(label_levels, dummy_levels)
  list(df = out_df, levels = out_levels)
}

# ----------------------------
# 8) Plot helper: faceted pies (identical facet cell size)
# ----------------------------
pie_panel <- function(df, sample_label_levels, panel_title, nrow_fixed, ncol_fixed) {
  
  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(panel_title, "\n(no samples)")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
    )
  }
  
  df <- df %>%
    mutate(Sample_label = factor(Sample_label, levels = sample_label_levels))
  
  ggplot(df, aes(x = 1, y = pct, fill = CellType4)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = pct_lab),
      position = position_stack(vjust = 0.5),
      size = 2.6,
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(values = celltype4_cols, drop = FALSE) +
    facet_wrap(
      ~ Sample_label,
      nrow = nrow_fixed,
      ncol = ncol_fixed,
      labeller = as_labeller(function(x) {
        ifelse(stringr::str_detect(x, "^\\.__dummy__"), "", x)
      })
    )+
    labs(title = panel_title, fill = "Cell Type") +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(size = 15, face = "bold", hjust = 0.5),
      strip.text      = element_text(size = 10.5, face = "bold", lineheight = 0.95),
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10),
      legend.key.size = unit(0.55, "cm"),
      panel.spacing   = unit(1.0, "lines")
    )
}


# ----------------------------
# 10) Build 3 panels WITH padding
# ----------------------------

# Pre-only
df_pre <- pie_df %>% filter(Patient %in% pre_only_patients, Timepoint == baseline_tp)
pad_pre  <- pad_to_grid(df_pre,  pre_only_order,  NROW_FIXED, NCOL_global)
p1_pre_only <- pie_panel(pad_pre$df, pad_pre$levels,
                         "Pre-treatment only", NROW_FIXED, NCOL_pre) +
  theme(strip.text = element_text(color = ifelse(grepl("^\\.__dummy__", as.character(pad_pre$levels)), NA, "black")))

# Paired
df_pair <- pie_df %>% filter(Patient %in% paired_patients,
                             Timepoint %in% c(baseline_tp, on_tps, post_tp))
pad_pair <- pad_to_grid(df_pair, paired_order,    NROW_FIXED, NCOL_global)
p2_paired <- pie_panel(pad_pair$df, pad_pair$levels,
                       "Paired pre-/on-treatment", NROW_FIXED, NCOL_pair)

# Post-only (this is the one that needed padding)
df_post <- pie_df %>% filter(Patient %in% post_only_patients, Timepoint == post_tp)
pad_post <- pad_to_grid(df_post, post_only_order, NROW_FIXED, NCOL_global)
p3_post_only <- pie_panel(pad_post$df, pad_post$levels,
                          "Post-treatment only", NROW_FIXED, NCOL_post) +
  theme(
    strip.text = element_text(
      # hide dummy facet titles
      colour = ifelse(grepl("^\\.__dummy__", as.character(pad_post$levels)), NA, "black")
    )
  )

# Combine: widths proportional to NCOL (facet cell size identical)
p_final <- (p1_pre_only + p2_paired + p3_post_only) +
  plot_layout(widths = c(1, 1, 1), guides = "collect")&
  theme(legend.position = "right")

print(p_final)




# ----------------------------
# 11) Export
# ----------------------------
facet_unit_in <- 1.9
legend_pad_in <- 3.0
fig_w <- facet_unit_in * (NCOL_pre + NCOL_pair + NCOL_post) + legend_pad_in
fig_h <- facet_unit_in * NROW_FIXED + 1.2

ggsave(filename = file.path(out_dir, "PieCellComposition__3cols__EqualPieSize.png"),
       p_final, width = fig_w, height = fig_h, units = "in", dpi = 900, bg = "white")

ggsave(filename = file.path(out_dir, "PieCellComposition__3cols__EqualPieSize.pdf"),
       p_final, width = fig_w, height = fig_h, units = "in", bg = "white")











# ----------------------------
# 9) Plot helper (square facet cells)
# ----------------------------
pie_panel <- function(df, sample_label_levels, panel_title, nrow_fixed, ncol_fixed) {
  
  if (nrow(df) == 0) {
    return(
      ggplot() + theme_void() +
        labs(title = paste0(panel_title, "\n(no samples)")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  }
  
  df <- df %>% mutate(Sample_label = factor(Sample_label, levels = sample_label_levels))
  
  ggplot(df, aes(x = 1, y = pct, fill = CellType4)) +
    geom_col(width = 1, color = "white", linewidth = 0.45) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = pct_lab),
      position = position_stack(vjust = 0.5),
      size = 2.6,
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(values = celltype4_cols, drop = FALSE) +
    facet_wrap(~ Sample_label, nrow = nrow_fixed, ncol = ncol_fixed) +
    labs(title = panel_title, fill = "Cell Type") +
    theme_void(base_size = 11) +
    theme(
      aspect.ratio = 1,  # <- forces each facet cell to be square
      plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
      strip.text   = element_text(size = 10.5, face = "bold", lineheight = 0.95),
      strip.background = element_blank(),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.55, "cm"),
      panel.spacing = unit(1.0, "lines")
    )
}

# ----------------------------
# 9) Build 3 panels (3 columns)
# ----------------------------
p1_pre_only <- pie_panel(
  df = pie_df %>% filter(Patient %in% pre_only_patients, Timepoint == baseline_tp),
  sample_label_levels = pre_only_order,
  panel_title = "Pre-treatment only",
  nrow_fixed = NROW_FIXED,
  ncol_fixed = NCOL_pre
)

p2_paired <- pie_panel(
  df = pie_df %>% filter(Patient %in% paired_patients,
                         Timepoint %in% c(baseline_tp, on_tps, post_tp)),
  sample_label_levels = paired_order,
  panel_title = "Paired pre-/on-treatment",
  nrow_fixed = NROW_FIXED,
  ncol_fixed = NCOL_pair
)

p3_post_only <- pie_panel(
  df = pie_df %>% filter(Patient %in% post_only_patients, Timepoint == post_tp),
  sample_label_levels = post_only_order,
  panel_title = "Post-treatment only",
  nrow_fixed = NROW_FIXED,
  ncol_fixed = NCOL_post
)

# Patchwork widths proportional to ncol -> identical pie sizes across panels
p_final <- (p1_pre_only + p2_paired + p3_post_only) +
  plot_layout(widths = c(NCOL_pre, NCOL_pair, NCOL_post), guides = "collect") &
  theme(legend.position = "right")

print(p_final)





# ----------------------------
# 10) Export: derive size from grid so pies remain consistent
# ----------------------------
facet_unit_in <- 1.9   # inches per facet cell (tweak if you want bigger pies)
legend_pad_in <- 3.0   # space for legend

fig_w <- facet_unit_in * (NCOL_pre + NCOL_pair + NCOL_post) + legend_pad_in
fig_h <- facet_unit_in * NROW_FIXED + 1.2

ggsave(filename = file.path(out_dir, "PieCellComposition__PreOnly__Paired__PostOnly__EqualPieSize.png"),
       p_final, width = fig_w, height = fig_h, units = "in", dpi = 900, bg = "white")

ggsave(filename = file.path(out_dir, "PieCellComposition__PreOnly__Paired__PostOnly__EqualPieSize.pdf"),
       p_final, width = fig_w, height = fig_h, units = "in", bg = "white")






# ----------------------------
# 8) Build the two panels and combine
# ----------------------------

p1_pre_only <- pie_panel(
  df = pie_df %>% filter(Patient %in% pre_only_patients, Timepoint == baseline_tp),
  sample_label_levels = pre_only_order,
  panel_title = "Pre-treatment only",
  ncol = 3
)

p2_paired <- pie_panel(
  df = pie_df %>% filter(Patient %in% paired_patients,
                         Timepoint %in% c(baseline_tp, on_tps, post_tp)),
  sample_label_levels = paired_order,
  panel_title = "Paired Pre + (OnT OR Post)",
  ncol = 4
)

p3_post_only <- pie_panel(
  df = pie_df %>% filter(Patient %in% post_only_patients,
                         Timepoint == post_tp),
  sample_label_levels = post_only_order,
  panel_title = "Post-treatment only",
  ncol = 1
)


# Combine into 3 columns with one shared legend
p_final <- (p1_pre_only + p2_paired + p3_post_only) +
  plot_layout(widths = c(1, 1, 1), guides = "collect") &
  theme(legend.position = "right")

print(p_final)


# ----------------------------
# 9) Export (high-res, journal-ready)
# ----------------------------
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig1def/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "PieCellComposition_BaselineOnly_vs_Paired.png"),
       p_final, width = 12, height = 6.5, units = "in", dpi = 900, bg = "white")

ggsave(filename = file.path(out_dir, "PieCellComposition_BaselineOnly_vs_Paired.pdf"),
       p_final, width = 12, height = 6.5, units = "in", bg = "white")
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig1def/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "PieCellComposition_BaselineOnly_vs_Paired_v2.png"),
       p_final, width = 12, height = 6.5, units = "in", dpi = 900, bg = "white")

ggsave(filename = file.path(out_dir, "PieCellComposition_BaselineOnly_vs_Paired_v2.pdf"),
       p_final, width = 12, height = 6.5, units = "in", bg = "white")





###############################################
## Pie charts per Sample_ID_renamed (3 columns) — aligned + equal pie size
## Columns:
##   1) Pre-treatment only
##   2) Paired pre-/on-treatment (includes Pre + OnT/C2T/C3T and/or Post)
##   3) Post-treatment only
##
## Key properties:
## - All three panels use the SAME facet grid: NROW_FIXED x NCOL_global
## - Panels are combined with equal widths (1:1:1) so pies stay identical size
## - Empty slots are padded with dummy facets BUT dummy labels are blanked via labeller
## - % labels: larger + adaptive (black on light slices, white on dark slices)
###############################################

suppressPackageStartupMessages({
})

# ------------------------------------------------------------
# 0) RM_ST assumed already loaded
# RM_ST <- readRDS(input_rds_1)
# ------------------------------------------------------------

# ----------------------------
# 0.1) Ensure Sample_ID_renamed exists
# ----------------------------
if (!"Sample_ID_renamed" %in% colnames(RM_ST@meta.data)) {
  message("Sample_ID_renamed not found. Creating it using Post->OnT logic...")
  
  keep_post <- c("LIB15_Post", "LIB17_Post")
  old_ids <- RM_ST$ID
  new_ids <- old_ids
  new_ids <- gsub("_Post$", "_OnT", new_ids)
  new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]
  RM_ST$Sample_ID_renamed <- new_ids
}

# ----------------------------
# 1) Create 4 broad categories
# ----------------------------
if (!"Spatial_CT_broad" %in% colnames(RM_ST@meta.data)) {
  stop("Spatial_CT_broad not found. Please create it first (as in your earlier script).")
}

RM_ST$CellType4 <- dplyr::case_when(
  RM_ST$Spatial_CT_broad %in% c("Malignant cells")    ~ "Malignant",
  RM_ST$Spatial_CT_broad %in% c("Fibroblasts")        ~ "Fibroblast",
  RM_ST$Spatial_CT_broad %in% c("Endothelial cells")  ~ "Endothelial",
  TRUE                                                ~ "Immune"
)

RM_ST$CellType4 <- factor(
  RM_ST$CellType4,
  levels = c("Endothelial", "Fibroblast", "Immune", "Malignant")
)

# ----------------------------
# 2) Parse patient + timepoint from Sample_ID_renamed
# ----------------------------
meta_df <- RM_ST@meta.data %>%
  mutate(
    Sample_ID_renamed = as.character(Sample_ID_renamed),
    Timepoint = str_extract(Sample_ID_renamed, "(?<=_)[^_]+$"),
    Patient   = str_replace(Sample_ID_renamed, "_[^_]+$", ""),
    Sample_label = paste0(Patient, "\n", Timepoint)
  )

# ----------------------------
# 2.5) Auto-detect baseline label ("BL" vs "Pre")
# ----------------------------
baseline_tp <- dplyr::case_when(
  "BL"  %in% meta_df$Timepoint ~ "BL",
  "Pre" %in% meta_df$Timepoint ~ "Pre",
  TRUE ~ NA_character_
)

if (is.na(baseline_tp)) {
  stop(
    "Could not find baseline timepoint label. Expected 'BL' or 'Pre'. Found: ",
    paste(sort(unique(meta_df$Timepoint)), collapse = ", ")
  )
}
message("Using baseline timepoint = ", baseline_tp)

# ----------------------------
# 2.6) Define on-treatment and post-treatment labels
# ----------------------------
on_candidates <- c("OnT", "C2T", "C3T")
on_tps <- intersect(on_candidates, unique(meta_df$Timepoint))
post_tp <- if ("Post" %in% unique(meta_df$Timepoint)) "Post" else NA_character_

message("Detected on-treatment timepoints: ", ifelse(length(on_tps)==0, "NONE", paste(on_tps, collapse=", ")))
message("Detected post-treatment timepoint: ", ifelse(is.na(post_tp), "NONE", post_tp))

# ----------------------------
# 3) Patient grouping into 3 columns (paired = Pre + (On or Post))
# ----------------------------
patient_tp <- meta_df %>%
  distinct(Patient, Timepoint) %>%
  group_by(Patient) %>%
  summarise(
    has_pre  = any(Timepoint == baseline_tp),
    has_on   = any(Timepoint %in% on_tps),
    has_post = ifelse(is.na(post_tp), FALSE, any(Timepoint == post_tp)),
    n_tp     = n(),
    .groups = "drop"
  )

pre_only_patients <- patient_tp %>%
  filter(has_pre, !has_on, !has_post, n_tp == 1) %>%
  pull(Patient)

paired_patients <- patient_tp %>%
  filter(has_pre, (has_on | has_post)) %>%
  pull(Patient)

post_only_patients <- patient_tp %>%
  filter(!has_pre, has_post, n_tp == 1) %>%
  pull(Patient)

# ----------------------------
# 4) Build per-sample proportions (+ visible labels)
# ----------------------------
pie_df <- meta_df %>%
  group_by(Patient, Timepoint, Sample_ID_renamed, Sample_label, CellType4) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample_ID_renamed) %>%
  mutate(
    pct = n / sum(n),
    pct_lab = ifelse(pct >= 0.05, percent(pct, accuracy = 1), "")
  ) %>%
  ungroup()

# ----------------------------
# 5) Ordering (Pre -> On -> Post)
# ----------------------------
tp_levels <- unique(c(baseline_tp, on_candidates, "Post"))
all_tp <- unique(pie_df$Timepoint)
tp_levels2 <- unique(c(tp_levels, setdiff(sort(all_tp), tp_levels)))

pie_df <- pie_df %>%
  mutate(
    Timepoint = factor(Timepoint, levels = tp_levels2),
    Patient   = factor(Patient, levels = sort(unique(Patient)))
  )

pre_only_order <- pie_df %>%
  filter(Patient %in% pre_only_patients, Timepoint == baseline_tp) %>%
  distinct(Sample_label, Patient, Timepoint) %>%
  arrange(Patient) %>%
  pull(Sample_label)

paired_order <- pie_df %>%
  filter(Patient %in% paired_patients,
         Timepoint %in% c(baseline_tp, on_tps, post_tp)) %>%
  distinct(Sample_label, Patient, Timepoint) %>%
  arrange(Patient, Timepoint) %>%
  pull(Sample_label)

post_only_order <- pie_df %>%
  filter(Patient %in% post_only_patients, Timepoint == post_tp) %>%
  distinct(Sample_label, Patient, Timepoint) %>%
  arrange(Patient) %>%
  pull(Sample_label)

# ----------------------------
# 6) Colors (journal-friendly)
# ----------------------------
celltype4_cols <- c(
  "Endothelial" = "#D9EF8B",
  "Fibroblast"  = "#C8BE9A",
  "Immune"      = "#2C6DB2",
  "Malignant"   = "#E67C73"
)
# ----------------------------
# 7) Fixed facet geometry (NO global grid, NO dummy padding)
# ----------------------------
NROW_FIXED <- 4

n_fac_pre  <- length(pre_only_order)
n_fac_pair <- length(paired_order)
n_fac_post <- length(post_only_order)

NCOL_pre  <- max(1, ceiling(n_fac_pre  / NROW_FIXED))
NCOL_pair <- max(1, ceiling(n_fac_pair / NROW_FIXED))
NCOL_post <- max(1, ceiling(n_fac_post / NROW_FIXED))

message("Facet grid (nrow fixed=", NROW_FIXED, "): ",
        "Pre-only ncol=", NCOL_pre, "; Paired ncol=", NCOL_pair, "; Post-only ncol=", NCOL_post)

# ----------------------------
# 8) Plot helper: faceted pies (no padding; equal pie sizes via widths + fixed nrow)
# ----------------------------
pie_panel <- function(df, sample_label_levels, panel_title, nrow_fixed, ncol_fixed) {
  
  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(panel_title, "\n(no samples)")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  }
  
  df <- df %>%
    mutate(Sample_label = factor(Sample_label, levels = sample_label_levels))
  
  ggplot(df, aes(x = 1, y = pct, fill = CellType4)) +
    geom_col(width = 1, color = "white", linewidth = 0.45) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = pct_lab),
      position = position_stack(vjust = 0.5),
      size = 2.6,
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(values = celltype4_cols, drop = FALSE) +
    facet_wrap(
      ~ Sample_label,
      nrow = nrow_fixed,
      ncol = ncol_fixed
    ) +
    labs(title = panel_title, fill = "Cell Type") +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(size = 15, face = "bold", hjust = 0.5),
      strip.text      = element_text(size = 10.5, face = "bold", lineheight = 0.95),
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10),
      legend.key.size = unit(0.55, "cm"),
      panel.spacing   = unit(1.0, "lines")
    )
}

# ----------------------------
# 9) Build the three panels (panel-specific ncol)
# ----------------------------
df_pre  <- pie_df %>% filter(Patient %in% pre_only_patients,  Timepoint == baseline_tp)
df_pair <- pie_df %>% filter(Patient %in% paired_patients,
                             Timepoint %in% c(baseline_tp, on_tps, post_tp))
df_post <- pie_df %>% filter(Patient %in% post_only_patients, Timepoint == post_tp)

p1_pre_only <- pie_panel(
  df = df_pre,
  sample_label_levels = pre_only_order,
  panel_title = "Pre-treatment only",
  nrow_fixed = NROW_FIXED,
  ncol_fixed = NCOL_pre
)

p2_paired <- pie_panel(
  df = df_pair,
  sample_label_levels = paired_order,
  panel_title = "Paired pre-/ (on/post treatment)",
  nrow_fixed = NROW_FIXED,
  ncol_fixed = NCOL_pair
)

p3_post_only <- pie_panel(
  df = df_post,
  sample_label_levels = post_only_order,
  panel_title = "Post-treatment only",
  nrow_fixed = NROW_FIXED,
  ncol_fixed = NCOL_post
)

# ----------------------------
# 10) Combine: widths proportional to panel ncol (THIS equalizes pie sizes)
# ----------------------------
p_final <- (p1_pre_only + p2_paired + p3_post_only) +
  plot_layout(widths = c(NCOL_pre, NCOL_pair, NCOL_post), guides = "collect") &
  theme(legend.position = "right")

print(p_final)

# ----------------------------
# 11) Export
# ----------------------------
facet_unit_in <- 1.9
legend_pad_in <- 3.0

fig_w <- facet_unit_in * (NCOL_pre + NCOL_pair + NCOL_post) + legend_pad_in
fig_h <- facet_unit_in * NROW_FIXED + 1.2

ggsave(filename = file.path(out_dir, "PieCellComposition__3cols__Compact_EqualPieSize.png"),
       p_final, width = fig_w, height = fig_h, units = "in", dpi = 900, bg = "white")

ggsave(filename = file.path(out_dir, "PieCellComposition__3cols__Compact_EqualPieSize.pdf"),
       p_final, width = fig_w, height = fig_h, units = "in", bg = "white")