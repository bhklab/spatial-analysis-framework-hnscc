# ============================================================
# Figure 5 — Additional analyses panels
# Publication polished script (paths parameterized, no setwd)
# ============================================================


# ---- shared helpers ----
# Assumes this script lives in the same folder as 00_common.R.
# If not, change the path below (e.g., source("code/04_part4/00_common.R")).
source("00_common.R")
theme_set(theme_pub())

suppressPackageStartupMessages({
  library(GSVA)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------
input_rds_1 <- "/part2/Cancer.rds"
input_rds_2 <- "/part2/Fibroblast_v2.RDS"
input_rds_3 <- "/part2/Macrophage.RDS"

# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "fig5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



#Figure 5
levels(factor(Cancer_v3$biopsy_location))
table(Cancer_v3$biopsy_location,useNA = "ifany")
########################################
## 0. Setup
########################################

# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig5/")  # removed (publication): use out_dir + file.path()
Cancer_v3 <- readRDS(input_rds_1)

## 1) Clean whitespace
Cancer_v3$biopsy_location_clean <- trimws(Cancer_v3$biopsy_location)

## 2) Map to lung / lymph_node / head_neck
Cancer_v3$biopsy_group <- dplyr::case_when(
  # LUNG
  Cancer_v3$biopsy_location_clean == "lung" ~ "lung",
  
  # LYMPH NODES
  Cancer_v3$biopsy_location_clean %in% c(
    "neck lymph node",
    "submandibular lymph node",
    "submental lymph node"
  ) ~ "lymph_node",
  
  # ALL OTHER HEAD/NECK SITES
  TRUE ~ "head_neck"
)

Cancer_v3$biopsy_group <- factor(
  Cancer_v3$biopsy_group,
  levels = c("lung", "lymph_node", "head_neck")
)

## QC
table(Cancer_v3$biopsy_location_clean)
table(Cancer_v3$biopsy_group, useNA = "ifany")


########################################
## 1. Sample ID renaming (Pre / OnT / Post)
########################################

keep_post <- c("LIB15_Post", "LIB17_Post")

old_ids <- Cancer_v3$ID
new_ids <- gsub("_Post$", "_OnT", old_ids)
new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]

Cancer_v3$Sample_ID_renamed <- factor(
  new_ids,
  levels = sort(unique(new_ids))
)
## Define responder sample IDs (by Sample_ID_renamed)
responder_ids <- c(
  "LIB04_Pre",
  "LIB17_Post",
  "LIB17_Pre",
  "S11_OnT",
  "S11_Pre",
  "S19_Pre",
  "S7_OnT",
  "S4_Pre"
)

## New column: Response_group
Cancer_v3$Response_group <- ifelse(
  Cancer_v3$Sample_ID_renamed %in% responder_ids,
  "Responder",
  "Non-responder"
)

Cancer_v3$Response_group <- factor(
  Cancer_v3$Response_group,
  levels = c("Responder", "Non-responder")
)

## QC
table(Cancer_v3$Sample_ID_renamed, Cancer_v3$Response_group)
table(Cancer_v3$Response_group)

levels(factor(Cancer_v3$Unsupervised_Clusters))
########################################
## 2. Malignant cluster renaming (0–4 -> biology labels)
########################################

cluster_name_map <- c(
  "Cluster 0" = "Cycling",
  "Cluster 1" = "TC",
  "Cluster 2" = "Neutrophil-inflamed",
  "Cluster 3" = "Fibrovascular",
  "Cluster 4" = "pEMT"
)

Cancer_v3$Malig_Lineage <- dplyr::recode(
  Cancer_v3$Unsupervised_Clusters,
  !!!cluster_name_map
)

Cancer_v3$Malig_Lineage <- factor(
  Cancer_v3$Malig_Lineage,
  levels = c("Cycling", "TC", "Neutrophil-inflamed", "Fibrovascular", "pEMT")
)

table(Cancer_v3$Unsupervised_Clusters, Cancer_v3$Malig_Lineage)

malig_cols <- c(
  "Cycling"             = "#D55E00",
  "TC"                  = "#0072B2",
  "Neutrophil-inflamed" = "#009E73",
  "Fibrovascular"       = "#CC79A7",
  "pEMT"                = "#F0E442"
)

########################################
## 3. Treatment phase (Pre / OnT / Post)
########################################

Cancer_v3$Treatment_phase_3 <- dplyr::case_when(
  grepl("_Pre$",  Cancer_v3$Sample_ID_renamed) ~ "Pre",
  grepl("_OnT$",  Cancer_v3$Sample_ID_renamed) ~ "OnT",
  grepl("_Post$", Cancer_v3$Sample_ID_renamed) ~ "Post",
  TRUE                                         ~ NA_character_
)

Cancer_v3$Treatment_phase_3 <- factor(
  Cancer_v3$Treatment_phase_3,
  levels = c("Pre", "OnT", "Post")
)

table(Cancer_v3$Sample_ID_renamed, Cancer_v3$Treatment_phase_3)
table(Cancer_v3$Treatment_phase_3, useNA = "ifany")

########################################
## 4. UMAP helpers and df_cancer (Pre/OnT only)
########################################

umap_cancer <- "umap.cca_Cancer_v3"  # check if needed

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

make_umap_df <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction)
  df  <- as.data.frame(emb)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df$cell <- rownames(emb)
  cbind(df, obj@meta.data[df$cell, , drop = FALSE])
}

df_cancer <- make_umap_df(Cancer_v3, umap_cancer) %>%
  dplyr::filter(Treatment_phase_3 %in% c("Pre", "OnT"))

df_cancer$Treatment_phase_3 <- droplevels(df_cancer$Treatment_phase_3)

xlim_cancer_pad <- pad_limits(range(df_cancer$UMAP_1), 0.03)
ylim_cancer_pad <- pad_limits(range(df_cancer$UMAP_2), 0.03)

########################################
## 5. Theme for UMAP
########################################

theme_umap_cancer <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      axis.line       = element_blank(),
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8),
      legend.position = "right",
      legend.key.size = unit(0.35, "cm"),
      plot.title      = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

########################################
## 6. Generic UMAP panel function
########################################

plot_umap_panel <- function(
    df,
    phase,            # "Pre" or "OnT"
    subset_idx,       # logical same length as nrow(df)
    pt.size   = 1.5,
    title_text = NULL
) {
  stopifnot(length(subset_idx) == nrow(df))
  
  df_sub <- df[subset_idx & df$Treatment_phase_3 == phase, , drop = FALSE]
  
  # Ensure Malig_Lineage exists and is factor
  if (!"Malig_Lineage" %in% colnames(df_sub)) {
    stop("Column 'Malig_Lineage' not found in df_sub.")
  }
  if (!is.factor(df_sub$Malig_Lineage)) {
    df_sub$Malig_Lineage <- factor(df_sub$Malig_Lineage,
                                   levels = levels(Cancer_v3$Malig_Lineage))
  }
  df_sub$Malig_Lineage <- droplevels(df_sub$Malig_Lineage)
  
  ggplot(df_sub, aes(x = UMAP_1, y = UMAP_2, color = Malig_Lineage)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_color_manual(values = malig_cols, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_cancer_pad,
      ylim   = ylim_cancer_pad,
      expand = FALSE
    ) +
    theme_umap_cancer() +
    labs(
      x     = NULL,
      y     = NULL,
      color = "Malignant lineage",
      title = title_text
    )
}

########################################
## 7. Indices for lung / TC / pEMT
########################################


table(df_cancer$id_patient, df_cancer$biopsy_group)

lung_idx <- df_cancer$biopsy_group == "lung"   # adjust value if needed
tc_idx   <- lung_idx & df_cancer$Malig_Lineage == "TC"
pemt_idx <- lung_idx & df_cancer$Malig_Lineage == "pEMT"

stopifnot(
  length(lung_idx) == nrow(df_cancer),
  length(tc_idx)   == nrow(df_cancer),
  length(pemt_idx) == nrow(df_cancer)
)

########################################
## 8. Build six panels (lung biopsies)
########################################

## Top row – all malignant lineages (lung)

p_lung_malig_pre <- plot_umap_panel(
  df         = df_cancer,
  phase      = "Pre",
  subset_idx = lung_idx,
  title_text = "Malignant clusters – Pre (Lung)"
)

p_lung_malig_on <- plot_umap_panel(
  df         = df_cancer,
  phase      = "OnT",
  subset_idx = lung_idx,
  title_text = "Malignant clusters – OnT (Lung)"
)

## Middle row – TC only, but same Malig_Lineage colors

p_lung_tc_pre <- plot_umap_panel(
  df         = df_cancer,
  phase      = "Pre",
  subset_idx = tc_idx,
  title_text = "TC cluster – Pre (Lung)"
)

p_lung_tc_on <- plot_umap_panel(
  df         = df_cancer,
  phase      = "OnT",
  subset_idx = tc_idx,
  title_text = "TC cluster – OnT (Lung)"
)

## Bottom row – pEMT only, same Malig_Lineage colors

p_lung_pemt_pre <- plot_umap_panel(
  df         = df_cancer,
  phase      = "Pre",
  subset_idx = pemt_idx,
  title_text = "pEMT cluster – Pre (Lung)"
)

p_lung_pemt_on <- plot_umap_panel(
  df         = df_cancer,
  phase      = "OnT",
  subset_idx = pemt_idx,
  title_text = "pEMT cluster – OnT (Lung)"
)

########################################
## 9. Combine 2×3 & save (single legend)
########################################

fig_lung <- (p_lung_malig_pre + p_lung_malig_on) /
  (p_lung_tc_pre    + p_lung_tc_on)    /
  (p_lung_pemt_pre  + p_lung_pemt_on) +
  plot_layout(guides = "collect")

fig_lung <- fig_lung & theme(legend.position = "right")

fig_lung

ggsave(
  filename = file.path(out_dir, "Fig_malignant_TC_pEMT_Lung_2x3.png"),
  plot     = fig_lung,
  width    = 8,
  height   = 7,
  units    = "in",
  dpi      = 900
)


########################################
## A. Indices for non-responders
########################################

# Quick sanity check
table(df_cancer$Response_group)

nonresp_idx      <- df_cancer$Response_group == "Non-responder"
nonresp_tc_idx   <- nonresp_idx & df_cancer$Malig_Lineage == "TC"
nonresp_pemt_idx <- nonresp_idx & df_cancer$Malig_Lineage == "pEMT"

stopifnot(
  length(nonresp_idx)      == nrow(df_cancer),
  length(nonresp_tc_idx)   == nrow(df_cancer),
  length(nonresp_pemt_idx) == nrow(df_cancer)
)



########################################
## B. Build six panels (non-responders)
########################################

## Top row – all malignant lineages within non-responders

p_nr_malig_pre <- plot_umap_panel(
  df         = df_cancer,
  phase      = "Pre",
  subset_idx = nonresp_idx,
  title_text = "Malignant clusters – Pre (Non-responders)"
)

p_nr_malig_on <- plot_umap_panel(
  df         = df_cancer,
  phase      = "OnT",
  subset_idx = nonresp_idx,
  title_text = "Malignant clusters – OnT (Non-responders)"
)

## Middle row – TC lineage only (still colored by Malig_Lineage)

p_nr_tc_pre <- plot_umap_panel(
  df         = df_cancer,
  phase      = "Pre",
  subset_idx = nonresp_tc_idx,
  title_text = "Tumore Core – Pre (Non-responders)"
)

p_nr_tc_on <- plot_umap_panel(
  df         = df_cancer,
  phase      = "OnT",
  subset_idx = nonresp_tc_idx,
  title_text = "Tumore Core – OnT (Non-responders)"
)

## Bottom row – pEMT lineage only

p_nr_pemt_pre <- plot_umap_panel(
  df         = df_cancer,
  phase      = "Pre",
  subset_idx = nonresp_pemt_idx,
  title_text = "pEMT – Pre (Non-responders)"
)

p_nr_pemt_on <- plot_umap_panel(
  df         = df_cancer,
  phase      = "OnT",
  subset_idx = nonresp_pemt_idx,
  title_text = "pEMT – OnT (Non-responders)"
)

########################################
## C. Combine 2×3 & save
########################################

fig_nonresp_panelA <- (p_nr_malig_pre + p_nr_malig_on) /
  (p_nr_tc_pre    + p_nr_tc_on)   +
  plot_layout(guides = "collect")

fig_nonresp_panelA <- fig_nonresp_panelA & theme(legend.position = "right")

fig_nonresp_panelA

fig_nonresp_panelB <- (p_nr_malig_pre + p_nr_malig_on) /
  (p_nr_pemt_pre  + p_nr_pemt_on) +
  plot_layout(guides = "collect")

fig_nonresp_panelB <- fig_nonresp_panelB & theme(legend.position = "right")

fig_nonresp_panelB
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig5/")  # removed (publication): use out_dir + file.path()
ggsave(
  filename = file.path(out_dir, "Fig_malignantfig_nonresp_panelA_2x2.png"),
  plot     = fig_nonresp_panelA,
  width    = 8,
  height   = 5,
  units    = "in",
  dpi      = 900
)

ggsave(
  filename = file.path(out_dir, "Fig_malignantfig_nonresp_panelB_2x2.png"),
  plot     = fig_nonresp_panelB,
  width    = 8,
  height   = 5,
  units    = "in",
  dpi      = 900
)


#Fibroblast


# Load fibroblast object
Fibroblast <- readRDS(input_rds_2)

########################################
## 1. Sample ID renaming (Pre / OnT / Post)
########################################

# IDs to preserve as true Post
keep_post <- c("LIB15_Post", "LIB17_Post")

old_ids <- Fibroblast$ID
new_ids <- gsub("_Post$", "_OnT", old_ids)
new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]

Fibroblast$Sample_ID_renamed <- factor(
  new_ids,
  levels = sort(unique(new_ids))
)

########################################
## 2. Treatment phase (Pre / OnT / Post)
########################################

Fibroblast$Treatment_phase_3 <- dplyr::case_when(
  grepl("_Pre$",  Fibroblast$Sample_ID_renamed) ~ "Pre",
  grepl("_OnT$",  Fibroblast$Sample_ID_renamed) ~ "OnT",
  grepl("_Post$", Fibroblast$Sample_ID_renamed) ~ "Post",
  TRUE                                          ~ NA_character_
)

Fibroblast$Treatment_phase_3 <- factor(
  Fibroblast$Treatment_phase_3,
  levels = c("Pre", "OnT", "Post")
)

table(Fibroblast$Sample_ID_renamed, Fibroblast$Treatment_phase_3)
table(Fibroblast$Treatment_phase_3, useNA = "ifany")

########################################
## 3. Responders vs Non-responders
########################################

responder_ids <- c(
  "LIB04_Pre",
  "LIB17_Post",
  "LIB17_Pre",
  "S11_OnT",
  "S11_Pre",
  "S19_Pre",
  "S7_OnT",
  "S7_Pre",
  "S4_Pre",
  "S4_OnT"
)

Fibroblast$Response_group <- ifelse(
  Fibroblast$Sample_ID_renamed %in% responder_ids,
  "Responder",
  "Non-responder"
)

Fibroblast$Response_group <- factor(
  Fibroblast$Response_group,
  levels = c("Responder", "Non-responder")
)

table(Fibroblast$Response_group)


########################################
## 4. Biopsy grouping
########################################

Fibroblast$biopsy_location_clean <- trimws(Fibroblast$biopsy_location)

Fibroblast$biopsy_group <- dplyr::case_when(
  # Lung
  Fibroblast$biopsy_location_clean == "lung" ~ "lung",
  
  # Lymph nodes
  Fibroblast$biopsy_location_clean %in% c(
    "neck lymph node",
    "submandibular lymph node",
    "submental lymph node"
  ) ~ "lymph_node",
  
  # Everything else = head/neck site
  TRUE ~ "head_neck"
)

Fibroblast$biopsy_group <- factor(
  Fibroblast$biopsy_group,
  levels = c("lung", "lymph_node", "head_neck")
)

table(Fibroblast$biopsy_location_clean)
table(Fibroblast$biopsy_group, useNA = "ifany")



########################################
## 5. Map Fibroblast clusters to biology labels
########################################

# If not already done, keep your existing C0–C5 mapping:
fibro_remap <- c(
  "Hypoxic/Invasion-associated CAFs (transitional)"        = "fibroblast_C0",
  "mCAFs"                                                  = "fibroblast_C1",
  "apCAFs"                                                 = "fibroblast_C2",
  "Hypoxic/Invasion-associated CAFs (high activation)"     = "fibroblast_C3",
  "myCAFs"                                                 = "fibroblast_C4",
  "Erythroid/Platelet-interacting CAFs"                    = "fibroblast_C5"
)

Fibroblast$Fibroblast_Cluster <- unname(
  fibro_remap[ Fibroblast$Fibroblast_SubCluster ]
)

Fibroblast$Fibroblast_Cluster <- factor(
  Fibroblast$Fibroblast_Cluster,
  levels = paste0("fibroblast_C", 0:5)
)

table(Fibroblast$Fibroblast_Cluster, useNA = "ifany")

# Now recode to biology-aware names
fibro_lineage_map <- c(
  "fibroblast_C0" = "non_activated_fibroblast",
  "fibroblast_C1" = "ecm_myCAF",
  "fibroblast_C2" = "IFNg_iCAF",
  "fibroblast_C3" = "IL_iCAF",
  "fibroblast_C4" = "acto_myCAF",
  "fibroblast_C5" = "erythroid_interacting"
)

Fibroblast$Fibro_Lineage <- dplyr::recode(
  Fibroblast$Fibroblast_Cluster,
  !!!fibro_lineage_map
)

Fibroblast$Fibro_Lineage <- factor(
  Fibroblast$Fibro_Lineage,
  levels = c(
    "non_activated_fibroblast",
    "ecm_myCAF",
    "IFNg_iCAF",
    "IL_iCAF",
    "acto_myCAF",
    "erythroid_interacting"
  )
)

table(Fibroblast$Fibroblast_Cluster, Fibroblast$Fibro_Lineage)


########################################
## 6. UMAP for Fibroblast – biology labels
########################################

Reductions(Fibroblast)
umap_fibro <- "umap.cca_Fibroblasts_v2"  # as you had

umap_fibro_df <- Embeddings(Fibroblast, umap_fibro)
xlim_fibro <- range(umap_fibro_df[, 1])
ylim_fibro <- range(umap_fibro_df[, 2])

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_fibro_pad <- pad_limits(xlim_fibro, 0.03)
ylim_fibro_pad <- pad_limits(ylim_fibro, 0.03)

# Colors for Fibro_Lineage (Okabe–Ito-ish, adjust as you like)
fibro_levels <- levels(Fibroblast$Fibro_Lineage)

fibro_cols <- setNames(
  c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#F0E442", "#56B4E9"),
  fibro_levels
)

theme_umap_fibro <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

fibro_umap <- DimPlot(
  Fibroblast,
  reduction = umap_fibro,
  group.by  = "Fibro_Lineage",
  pt.size   = 0.8,
  raster    = FALSE
) +
  scale_color_manual(values = fibro_cols, drop = FALSE) +
  coord_cartesian(
    xlim   = xlim_fibro_pad,
    ylim   = ylim_fibro_pad,
    expand = FALSE
  ) +
  theme_umap_fibro() +
  labs(
    x     = NULL,
    y     = NULL,
    color = "Fibroblast lineage",
    title = "Fibroblast UMAP (non_activated / myCAF / iCAF / acto_myCAF / erythroid_interacting)"
  )

fibro_umap


########################################
## 4. UMAP helpers & df_fibro (Pre/OnT only)
########################################

umap_fibro <- "umap.cca_Fibroblasts_v2"  # as before

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

make_umap_df <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction)
  df  <- as.data.frame(emb)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df$cell <- rownames(emb)
  cbind(df, obj@meta.data[df$cell, , drop = FALSE])
}

df_fibro <- make_umap_df(Fibroblast, umap_fibro) %>%
  dplyr::filter(Treatment_phase_3 %in% c("Pre", "OnT"))

df_fibro$Treatment_phase_3 <- droplevels(df_fibro$Treatment_phase_3)

xlim_fibro_pad <- pad_limits(range(df_fibro$UMAP_1), 0.03)
ylim_fibro_pad <- pad_limits(range(df_fibro$UMAP_2), 0.03)

########################################
## 5. Colors & theme
########################################

fibro_levels <- levels(Fibroblast$Fibro_Lineage)

fibro_cols <- setNames(
  c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#F0E442", "#56B4E9"),
  fibro_levels
)

theme_umap_fibro <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

########################################
## 6. Fibro UMAP panel function
########################################

plot_fibro_panel <- function(
    df,
    phase,
    subset_idx,
    pt.size   = 1.5,
    title_text = NULL,
    drop_ery = TRUE    # NEW ARGUMENT
) {
  stopifnot(length(subset_idx) == nrow(df))
  
  df_sub <- df[subset_idx & df$Treatment_phase_3 == phase, , drop = FALSE]
  
  if (drop_ery) {
    df_sub <- df_sub[df_sub$Fibro_Lineage != "erythroid_interacting", , drop = FALSE]
  }
  
  df_sub$Fibro_Lineage <- droplevels(df_sub$Fibro_Lineage)
  
  ggplot(df_sub, aes(x = UMAP_1, y = UMAP_2, color = Fibro_Lineage)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_color_manual(values = fibro_cols, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_fibro_pad,
      ylim   = ylim_fibro_pad,
      expand = FALSE
    ) +
    theme_umap_fibro() +
    labs(
      x     = NULL,
      y     = NULL,
      color = "Fibroblast lineage",
      title = title_text
    )
}

########################################
## 7. Indices for responders / acto_myCAF
########################################

resp_idx      <- df_fibro$Response_group == "Responder"
acto_idx      <- resp_idx & df_fibro$Fibro_Lineage == "acto_myCAF"

# For bottom row, we can show acte_myCAF in lung responders as an example
acto_lung_idx <- resp_idx &
  df_fibro$Fibro_Lineage == "acto_myCAF" &
  df_fibro$biopsy_group == "lung"

stopifnot(
  length(resp_idx)      == nrow(df_fibro),
  length(acto_idx)      == nrow(df_fibro),
  length(acto_lung_idx) == nrow(df_fibro)
)

########################################
## 8. Build six panels – RESPONDERS
########################################

## Top row – all fibroblast lineages in responders

p_resp_fibro_pre <- plot_fibro_panel(
  df         = df_fibro,
  phase      = "Pre",
  subset_idx = resp_idx,
  title_text = "CAFs – Pre (Responders)",
  drop_ery   = TRUE        # <--- exclude erythroid_interacting
  
)

p_resp_fibro_on <- plot_fibro_panel(
  df         = df_fibro,
  phase      = "OnT",
  subset_idx = resp_idx,
  title_text = "CAFs – OnT (Responders)",
  drop_ery   = TRUE        # <--- exclude erythroid_interacting
  
)

## Middle row – acto_myCAF in responders (all sites)

p_resp_acto_pre <- plot_fibro_panel(
  df         = df_fibro,
  phase      = "Pre",
  subset_idx = acto_idx,
  title_text = "acto_myCAF – Pre (Responders)",
  drop_ery   = TRUE        # <--- exclude erythroid_interacting
  
)

p_resp_acto_on <- plot_fibro_panel(
  df         = df_fibro,
  phase      = "OnT",
  subset_idx = acto_idx,
  title_text = "acto_myCAF – OnT (Responders)",
  drop_ery   = TRUE        # <--- exclude erythroid_interacting
  
)

########################################
## 9. Combine 2×3 & save
########################################

fig_fibro_resp <- (p_resp_fibro_pre      + p_resp_fibro_on)      /
  (p_resp_acto_pre       + p_resp_acto_on)       +
  plot_layout(guides = "collect")

fig_fibro_resp <- fig_fibro_resp & theme(legend.position = "right")

fig_fibro_resp
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig5/")  # removed (publication): use out_dir + file.path()
ggsave(
  filename = file.path(out_dir, "Fig_Fibroblast_acto_myCAF_Responders_2x2.png"),
  plot     = fig_fibro_resp,
  width    = 8,
  height   = 5,
  units    = "in",
  dpi      = 900
)






#Macrophages
############################################################
## 0. Setup & load
############################################################


Macrophage <- readRDS(input_rds_3)

############################################################
## 1. Sample ID renaming (Pre / OnT / Post)
############################################################

keep_post <- c("LIB15_Post", "LIB17_Post")

old_ids <- Macrophage$ID
new_ids <- gsub("_Post$", "_OnT", old_ids)
new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]

Macrophage$Sample_ID_renamed <- factor(
  new_ids,
  levels = sort(unique(new_ids))
)

############################################################
## 2. Treatment phase (Pre / OnT / Post)
############################################################

Macrophage$Treatment_phase_3 <- dplyr::case_when(
  grepl("_Pre$",  Macrophage$Sample_ID_renamed) ~ "Pre",
  grepl("_OnT$",  Macrophage$Sample_ID_renamed) ~ "OnT",
  grepl("_Post$", Macrophage$Sample_ID_renamed) ~ "Post",
  TRUE                                          ~ NA_character_
)

Macrophage$Treatment_phase_3 <- factor(
  Macrophage$Treatment_phase_3,
  levels = c("Pre", "OnT", "Post")
)

table(Macrophage$Sample_ID_renamed, Macrophage$Treatment_phase_3)
table(Macrophage$Treatment_phase_3, useNA = "ifany")

############################################################
## 3. Responders vs Non-responders (same definition as before)
############################################################

responder_ids <- unique(c(
  "LIB04_Pre",
  "LIB17_Post",
  "LIB17_Pre",
  "S11_OnT",
  "S11_Pre",
  "S19_Pre",
  "S7_Pre",
  "S7_OnT",
  "S4_Pre",
  "S4_OnT"
))

Macrophage$Response_group <- ifelse(
  Macrophage$Sample_ID_renamed %in% responder_ids,
  "Responder",
  "Non-responder"
)

Macrophage$Response_group <- factor(
  Macrophage$Response_group,
  levels = c("Responder", "Non-responder")
)

table(Macrophage$Response_group)

############################################################
## 4. Biopsy grouping (lung / lymph_node / head_neck)
############################################################

Macrophage$biopsy_location_clean <- trimws(Macrophage$biopsy_location)

Macrophage$biopsy_group <- dplyr::case_when(
  Macrophage$biopsy_location_clean == "lung" ~ "lung",
  Macrophage$biopsy_location_clean %in% c(
    "neck lymph node",
    "submandibular lymph node",
    "submental lymph node"
  ) ~ "lymph_node",
  TRUE ~ "head_neck"
)

Macrophage$biopsy_group <- factor(
  Macrophage$biopsy_group,
  levels = c("lung", "lymph_node", "head_neck")
)

table(Macrophage$biopsy_location_clean)
table(Macrophage$biopsy_group, useNA = "ifany")

############################################################
## 5. HPV status collapse (for HPV– figure later)
############################################################

if ("hpv_status" %in% colnames(Macrophage@meta.data)) {
  Macrophage$HPV_status_collapsed <- as.character(Macrophage$hpv_status)
  Macrophage$HPV_status_collapsed[
    Macrophage$HPV_status_collapsed == "Non-tested"
  ] <- "Negative"
  
  Macrophage$HPV_status_collapsed <- factor(
    Macrophage$HPV_status_collapsed,
    levels = c("Negative", "Positive")
  )
  
  table(Macrophage$HPV_status_collapsed, useNA = "ifany")
}

############################################################
## 6. Map Macrophage clusters to C0–C2 and biology labels
############################################################

table(Macrophage$Macrophage_Identity, useNA = "ifany")

macro_remap <- c(
  "M2-like/ECM-remodeling TAMs"   = "macrophage_C0",
  "M1-like/Inflammatory TAMs"     = "macrophage_C1",
  "tTAMs"                         = "macrophage_C2"
)

Macrophage$Macrophage_Cluster <- unname(
  macro_remap[ Macrophage$Macrophage_Identity ]
)

Macrophage$Macrophage_Cluster <- factor(
  Macrophage$Macrophage_Cluster,
  levels = paste0("macrophage_C", 0:2)
)

table(Macrophage$Macrophage_Cluster, useNA = "ifany")

# Biology-aware lineages (short names for plotting)
macro_lineage_map <- c(
  "macrophage_C0" = "IFN_TAM",
  "macrophage_C1" = "inflammatory_TAM",
  "macrophage_C2" = "monocyte"
)


Macrophage$Macro_Lineage <- dplyr::recode(
  Macrophage$Macrophage_Cluster,
  !!!macro_lineage_map
)

Macrophage$Macro_Lineage <- factor(
  Macrophage$Macro_Lineage,
  levels = c("IFN_TAM", "inflammatory_TAM", "monocyte")
)

table(Macrophage$Macrophage_Cluster, Macrophage$Macro_Lineage)

############################################################
## 7. UMAP for Macrophage – lineage labels
############################################################

Reductions(Macrophage)
umap_macro <- "umap.cca_macro"

umap_macro_df <- Embeddings(Macrophage, umap_macro)
xlim_macro <- range(umap_macro_df[, 1])
ylim_macro <- range(umap_macro_df[, 2])

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_macro_pad <- pad_limits(xlim_macro, 0.03)
ylim_macro_pad <- pad_limits(ylim_macro, 0.03)

macro_levels <- levels(Macrophage$Macro_Lineage)

macro_cols <- c(
  "IFN_TAM"          = "#CC6677",  # reddish
  "inflammatory_TAM" = "#882255",  # deep wine
  "monocyte"                = "#44AA99"   # teal
)

theme_umap_macro <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

macro_umap <- DimPlot(
  Macrophage,
  reduction = umap_macro,
  group.by  = "Macro_Lineage",
  pt.size   = 0.8,
  raster    = FALSE
) +
  scale_color_manual(
    values = macro_cols,
    drop   = FALSE
  ) +
  coord_cartesian(
    xlim   = xlim_macro_pad,
    ylim   = ylim_macro_pad,
    expand = FALSE
  ) +
  theme_umap_macro() +
  labs(
    x     = NULL,
    y     = NULL,
    color = "Macrophage lineage",
    title = "Macrophage UMAP (M2_ECM_TAM / M1_inflammatory_TAM / tTAM)"
  )

macro_umap

############################################################
## 8. UMAP helpers & df_macro (Pre/OnT only)
############################################################

make_umap_df <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction)
  df  <- as.data.frame(emb)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df$cell <- rownames(emb)
  cbind(df, obj@meta.data[df$cell, , drop = FALSE])
}

df_macro <- make_umap_df(Macrophage, umap_macro) %>%
  dplyr::filter(Treatment_phase_3 %in% c("Pre", "OnT"))

df_macro$Treatment_phase_3 <- droplevels(df_macro$Treatment_phase_3)

xlim_macro_pad <- pad_limits(range(df_macro$UMAP_1), 0.03)
ylim_macro_pad <- pad_limits(range(df_macro$UMAP_2), 0.03)

############################################################
## 9. Generic Macrophage UMAP panel function
############################################################

theme_umap_macro <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

plot_macro_panel <- function(
    df,
    phase,            # "Pre" or "OnT"
    subset_idx,       # logical same length as nrow(df)
    pt.size   = 1.2,
    title_text = NULL
) {
  stopifnot(length(subset_idx) == nrow(df))
  
  df_sub <- df[subset_idx & df$Treatment_phase_3 == phase, , drop = FALSE]
  
  if (!"Macro_Lineage" %in% colnames(df_sub)) {
    stop("Column 'Macro_Lineage' not found in df_sub.")
  }
  if (!is.factor(df_sub$Macro_Lineage)) {
    df_sub$Macro_Lineage <- factor(df_sub$Macro_Lineage,
                                   levels = levels(Macrophage$Macro_Lineage))
  }
  df_sub$Macro_Lineage <- droplevels(df_sub$Macro_Lineage)
  
  ggplot(df_sub, aes(x = UMAP_1, y = UMAP_2, color = Macro_Lineage)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_color_manual(values = macro_cols, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_macro_pad,
      ylim   = ylim_macro_pad,
      expand = FALSE
    ) +
    theme_umap_macro() +
    labs(
      x     = NULL,
      y     = NULL,
      color = "Macrophage lineage",
      title = title_text
    )
}




## If HPV_status_collapsed not made yet, create it
if (!"HPV_status_collapsed" %in% colnames(Macrophage@meta.data) &&
    "hpv_status" %in% colnames(Macrophage@meta.data)) {
  Macrophage$HPV_status_collapsed <- as.character(Macrophage$hpv_status)
  Macrophage$HPV_status_collapsed[
    Macrophage$HPV_status_collapsed == "Non-tested"
  ] <- "Negative"
  
  Macrophage$HPV_status_collapsed <- factor(
    Macrophage$HPV_status_collapsed,
    levels = c("Negative", "Positive")
  )
}

table(Macrophage$HPV_status_collapsed, useNA = "ifany")

## Helper to build df_macro with UMAP coords + metadata
umap_macro <- "umap.cca_macro"

make_umap_df <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction)
  df  <- as.data.frame(emb)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df$cell <- rownames(emb)
  cbind(df, obj@meta.data[df$cell, , drop = FALSE])
}

df_macro <- make_umap_df(Macrophage, umap_macro) %>%
  dplyr::filter(Treatment_phase_3 %in% c("Pre", "OnT"))

df_macro$Treatment_phase_3 <- droplevels(df_macro$Treatment_phase_3)

## UMAP limits with padding
pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_macro_pad <- pad_limits(range(df_macro$UMAP_1), 0.03)
ylim_macro_pad <- pad_limits(range(df_macro$UMAP_2), 0.03)


## Macro lineage colors (as before)
macro_cols <- c(
  "IFN_TAM"          = "#CC6677",  # reddish
  "inflammatory_TAM" = "#882255",  # deep wine
  "monocyte"                = "#44AA99"   # teal
)

theme_umap_macro <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

## 3A. Macrophage panel (colors = Macro_Lineage)

plot_macro_panel <- function(
    df,
    phase,            # "Pre" or "OnT"
    subset_idx,       # logical same length as nrow(df)
    pt.size   = 1.2,
    title_text = NULL
) {
  stopifnot(length(subset_idx) == nrow(df))
  
  df_sub <- df[subset_idx & df$Treatment_phase_3 == phase, , drop = FALSE]
  
  if (!"Macro_Lineage" %in% colnames(df_sub)) {
    stop("Column 'Macro_Lineage' not found in df_sub.")
  }
  if (!is.factor(df_sub$Macro_Lineage)) {
    df_sub$Macro_Lineage <- factor(df_sub$Macro_Lineage,
                                   levels = levels(Macrophage$Macro_Lineage))
  }
  df_sub$Macro_Lineage <- droplevels(df_sub$Macro_Lineage)
  
  ggplot(df_sub, aes(x = UMAP_1, y = UMAP_2, color = Macro_Lineage)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_color_manual(values = macro_cols, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_macro_pad,
      ylim   = ylim_macro_pad,
      expand = FALSE
    ) +
    theme_umap_macro() +
    labs(
      x     = NULL,
      y     = NULL,
      color = "Macrophage lineage",
      title = title_text
    )
}

## 3B. Monocyte panel (colors = monocyte subset column)
##  >>>>> CHANGE THIS to your real column name if different
monocyte_col <- "monocyte"

plot_monocyte_panel <- function(
    df,
    phase,
    subset_idx,
    pt.size   = 1.2,
    title_text = NULL
) {
  stopifnot(length(subset_idx) == nrow(df))
  
  df_sub <- df[subset_idx & df$Treatment_phase_3 == phase, , drop = FALSE]
  
  if (!monocyte_col %in% colnames(df_sub)) {
    stop(sprintf("Column '%s' not found in df_sub.", monocyte_col))
  }
  
  # Keep only cells with non-missing monocyte annotation
  df_sub <- df_sub[!is.na(df_sub[[monocyte_col]]), , drop = FALSE]
  
  df_sub[[monocyte_col]] <- droplevels(factor(df_sub[[monocyte_col]]))
  mono_levels <- levels(df_sub[[monocyte_col]])
  
  mono_cols <- setNames(
    grDevices::hcl.colors(length(mono_levels), palette = "Dark3"),
    mono_levels
  )
  
  ggplot(df_sub, aes(x = UMAP_1, y = UMAP_2, color = .data[[monocyte_col]])) +
    geom_point(size = pt.size, stroke = 0) +
    scale_color_manual(values = mono_cols, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_macro_pad,
      ylim   = ylim_macro_pad,
      expand = FALSE
    ) +
    theme_umap_macro() +
    labs(
      x     = NULL,
      y     = NULL,
      color = "Monocyte subset",
      title = title_text
    )
}


## HPV-negative indices
hpv_neg_idx <- df_macro$HPV_status_collapsed == "Negative"

########################################
## Colors & theme (once)
########################################

macro_cols <- c(
  "IFN_TAM"          = "#CC6677",  # reddish
  "inflammatory_TAM" = "#882255",  # deep wine
  "monocyte"         = "#44AA99"   # teal
)

theme_umap_macro <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

plot_macro_panel <- function(
    df,
    phase,            # "Pre" or "OnT"
    subset_idx,       # logical same length as nrow(df)
    pt.size   = 1.2,
    title_text = NULL
) {
  stopifnot(length(subset_idx) == nrow(df))
  
  df_sub <- df[subset_idx & df$Treatment_phase_3 == phase, , drop = FALSE]
  
  if (!"Macro_Lineage" %in% colnames(df_sub)) {
    stop("Column 'Macro_Lineage' not found in df_sub.")
  }
  if (!is.factor(df_sub$Macro_Lineage)) {
    df_sub$Macro_Lineage <- factor(
      df_sub$Macro_Lineage,
      levels = levels(Macrophage$Macro_Lineage)
    )
  }
  df_sub$Macro_Lineage <- droplevels(df_sub$Macro_Lineage)
  
  ggplot(df_sub, aes(x = UMAP_1, y = UMAP_2, color = Macro_Lineage)) +
    geom_point(size = pt.size, stroke = 0) +
    scale_color_manual(values = macro_cols, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_macro_pad,
      ylim   = ylim_macro_pad,
      expand = FALSE
    ) +
    theme_umap_macro() +
    labs(
      x     = NULL,
      y     = NULL,
      color = "Macrophage lineage",
      title = title_text
    )
}

########################################
## HPV-negative indices
########################################

hpv_neg_idx <- df_macro$HPV_status_collapsed == "Negative"
stopifnot(length(hpv_neg_idx) == nrow(df_macro))

# Top row: all macrophage lineages, HPV–
macro_hpvneg_idx <- hpv_neg_idx

# Bottom row: monocyte lineage only, HPV–
mono_hpvneg_idx <- hpv_neg_idx & df_macro$Macro_Lineage == "monocyte"
stopifnot(length(mono_hpvneg_idx) == nrow(df_macro))

########################################
## Build four panels – HPV-negative
########################################

## Top row – all TAMs (all lineages), HPV–
p_hpvneg_macro_pre <- plot_macro_panel(
  df         = df_macro,
  phase      = "Pre",
  subset_idx = macro_hpvneg_idx,
  title_text = "TAMs – Pre (HPV–)"
)

p_hpvneg_macro_on <- plot_macro_panel(
  df         = df_macro,
  phase      = "OnT",
  subset_idx = macro_hpvneg_idx,
  title_text = "TAMs – OnT (HPV–)"
)

## Bottom row – monocyte lineage only, HPV–
p_hpvneg_mono_pre <- plot_macro_panel(
  df         = df_macro,
  phase      = "Pre",
  subset_idx = mono_hpvneg_idx,
  title_text = "Monocyte – Pre (HPV–)"
)

p_hpvneg_mono_on <- plot_macro_panel(
  df         = df_macro,
  phase      = "OnT",
  subset_idx = mono_hpvneg_idx,
  title_text = "Monocyte – OnT (HPV–)"
)

########################################
## Combine 2×2 & save
########################################

fig_macro_hpvneg <- (p_hpvneg_macro_pre + p_hpvneg_macro_on) /
  (p_hpvneg_mono_pre  + p_hpvneg_mono_on)  +
  plot_layout(guides = "collect")

fig_macro_hpvneg <- fig_macro_hpvneg & theme(legend.position = "right")

fig_macro_hpvneg
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig5/")  # removed (publication): use out_dir + file.path()
ggsave(
  filename = file.path(out_dir, "Fig_Macrophage_Monocyte_HPVneg_2x2.png"),
  plot     = fig_macro_hpvneg,
  width    = 8,
  height   = 5,
  units    = "in",
  dpi      = 900
)