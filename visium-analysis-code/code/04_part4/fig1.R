# ============================================================
# Figure 1 — Cohort overview (spatial + UMAP)
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
  library(stringr)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------
input_rds_1 <- "/part2/RM_ST.rds"

# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "fig1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################
## 0. Setup environment
###############################################

# Load libraries

# Load object
RM_ST <- readRDS(input_rds_1)



#Figure 1def
# Where you want to save figures
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/")  # removed (publication): use out_dir + file.path()
set.seed(1234)  # for reproducibility if you ever rerun UMAP

# Keep clustering neighbors/graph as-is (dims 1:15), just make a "visual UMAP"
RM_ST <- RunUMAP(
  RM_ST,
  reduction     = "integrated.dr",
  dims          = 1:5,                  # <- fewer PCs for viz
  reduction.name = "umap.integrated.vis"
)
ElbowPlot(RM_ST)
p_vis <- DimPlot(
  RM_ST,
  reduction = "umap.integrated.vis",
  group.by  = "Pathologist_Annotation",
  raster    = F,
  pt.size   = 1
)

###############################################
## 1. Quick sanity checks
###############################################

levels(factor(RM_ST$Spatial_CT))
table(RM_ST$Spatial_CT, useNA = "ifany")

levels(factor(RM_ST$Pathologist_Annotation))
table(RM_ST$Pathologist_Annotation, useNA = "ifany")

levels(factor(RM_ST$ID))
table(RM_ST$ID, useNA = "ifany")

###############################################
## 2. Clean + order pathologist annotations
###############################################

# Make sure levels are consistent and ordered for plotting
RM_ST$Pathologist_Annotation <- factor(
  RM_ST$Pathologist_Annotation,
  levels = c("SCC","Non-SCC", "Out-of-slide")
)

# Color palette (publication-friendly)
ann_cols <- c(
  "Non-SCC"      = "#08519C",   # green
  "SCC"          = "#C51B7D",   # magenta-red (strong tumor color)
  "Out-of-slide" = "#999999"    # dark neutral gray
)


mal_cols <- c(
  "Malignant spots"     = "#FF0000",  # deep magenta (tumor)
  "Non-Malignant spots" = "darkgrey"   # green (stromal/immune)
)


#We need some help with pin-pointing the 9 OnT samples and 2 Post-Treatment Samles
#LIB15_Post => LIB15_OnT
#LIB17_Post => LIB17_OnT

DimPlot(RM_ST,reduction = "umap.integrated", group.by = "Pathologist_Annotation")
DimPlot(RM_ST,reduction = "umap.integrated", group.by = "Pathologist_Annotation")




###############################################
## 1. Define malignant Spatial_CT categories
###############################################

malignant_clusters <- c(
  "TC_Epithelial_Differentiated",
  "Neuroendocrine/Stem-like Niche",
  "LE_ECM_Remodeling",
  "LE_Pro_Inflammatory",
  "Hypoxic/Angiogenic Niche"
)

RM_ST$Consensus_Malignancy <- ifelse(
  RM_ST$Spatial_CT %in% malignant_clusters,
  "Malignant spots",
  "Non-Malignant spots"
)

# Turn into an ordered factor
RM_ST$Consensus_Malignancy <- factor(
  RM_ST$Consensus_Malignancy,
  levels = c("Malignant spots", "Non-Malignant spots")
)

table(RM_ST$Consensus_Malignancy)


############################################################
## 0. UMAP limits (run once after UMAP is computed)
############################################################

umap_df <- Embeddings(RM_ST, "umap.integrated.vis")

xlim_fixed <- range(umap_df[, 1])
ylim_fixed <- range(umap_df[, 2])

xlim_fixed
ylim_fixed

############################################################
## 1. Color palettes
############################################################

# Pathologist annotation colors
ann_cols <- c(
  "Non-SCC"      = "darkblue",   # or e.g. "#08519C" if you want hex
  "SCC"          = "#C51B7D",    # magenta-red tumor
  "Out-of-slide" = "#999999"     # grey
)

# Consensus malignancy colors with brighter red
mal_cols <- c(
  "Malignant spots"     = "#D41114",   # bright strong red
  "Non-Malignant spots" = "#999999"    # green
)

############################################################
## 2. Theme for PNG figures (no box, legend ON)
############################################################

theme_umap <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),   # remove axis titles
      axis.text        = element_blank(),   # remove axis tick labels
      axis.ticks       = element_blank(),   # remove tick marks
      axis.line        = element_blank(),   # remove axis lines
      
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      
      plot.title       = element_text(size = 11, face = "bold", 
                                      hjust = 0.5)
    )
}


############################################################
## 3. Standardized UMAP function for PNGs (with legend)
############################################################

umap_png <- function(group_by, colors, title_text) {
  
  DimPlot(
    RM_ST,
    reduction = "umap.integrated.vis",
    group.by  = group_by,
    pt.size   = 0.8,     # 1.5 is very chunky; 0.4–0.8 is usually ideal
    raster    = F
  ) +
    scale_color_manual(values = colors, drop = FALSE) +
    
    # Fixed coordinate range → identical panel framing
    coord_cartesian(
      xlim   = xlim_fixed,
      ylim   = ylim_fixed,
      expand = FALSE
    ) +
    
    theme_umap() +
    
    labs(
      x     = "UMAP 1",
      y     = "UMAP 2",
      color = title_text,
      title = title_text
    )
}

############################################################
## 4. Build the UMAPs
############################################################

p_path_png <- umap_png(
  group_by   = "Pathologist_Annotation",
  colors     = ann_cols,
  title_text = "Pathologist Annotation"
)
print(p_path_png)

p_consensus_png <- umap_png(
  group_by   = "Consensus_Malignancy",
  colors     = mal_cols,
  title_text = "Consensus-Based Malignancy"
)
print(p_consensus_png)

############################################################
## 5. High-resolution PNG export (identical size)
############################################################

ggsave(filename = file.path(out_dir, "UMAP_PathologistAnnotation.png"),
  p_path_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "UMAP_ConsensusMalignancy.png"),
  p_consensus_png,
  width  = 7.5,
  height = 5,
  units  = "in",
  dpi    = 900
)

p1 <- SpatialDimPlot(
  RM_ST,
  group.by = "Consensus_Malignancy",
  images = "spatial.6",
  cols = mal_cols,
  pt.size.factor = 5,
  combine = FALSE
)[[1]]

p2 <- SpatialDimPlot(
  RM_ST,
  group.by = "Pathologist_Annotation",
  images = "spatial.6",
  cols = ann_cols,
  pt.size.factor = 5,
  combine = FALSE
)[[1]]


p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")


ggsave(filename = file.path(out_dir, "Spatial_spatial6_ConsensusMalignancy_noLegend.png"),
  p1,
  width = 7,
  height = 5,
  units = "in",
  dpi = 900
)

ggsave(filename = file.path(out_dir, "Spatial_spatial6_PathologistAnnotation_noLegend.png"),
  p2,
  width = 7,
  height = 5,
  units = "in",
  dpi = 900
)


#Figures
malignant_clusters <- c(
  "TC_Epithelial_Differentiated",
  "Neuroendocrine/Stem-like Niche",
  "LE_ECM_Remodeling",
  "LE_Pro_Inflammatory",
  "Hypoxic/Angiogenic Niche"
)

caf_clusters <- c(
  "apCAFs",
  "mCAFs",
  "myCAFs",
  "Hypoxic/Invasion-associated CAFs (high activation)",
  "Hypoxic/Invasion-associated CAFs (transitional)",
  "Erythroid/Platelet-interacting CAFs"
)

tam_clusters <- c(
  "M1-like/Inflammatory TAMs",
  "M2-like/ECM-remodeling TAMs",
  "tTAMs"
)

bplasma_clusters <- c(
  "Mature B cell/Plasma cell",
  "Mucosal Plasma B cells",
  "tBcells"
)


RM_ST$Spatial_CT_broad <- case_when(
  RM_ST$Spatial_CT %in% malignant_clusters ~ "Malignant cells",
  RM_ST$Spatial_CT %in% caf_clusters       ~ "Fibroblasts",
  RM_ST$Spatial_CT %in% tam_clusters       ~ "Macrophages",
  RM_ST$Spatial_CT %in% bplasma_clusters   ~ "B cells",
  
  # everything else retains original identity
  TRUE ~ RM_ST$Spatial_CT
)

final_levels <- c(
  "Malignant cells",
  "Fibroblasts",
  "Macrophages",
  "B cells",
  "T cells",
  "Dendritic cells",
  "Endothelial cells",
  "Terminal keratinocytes",
  "Muscle cells",
  "Salivary gland cells"
)

RM_ST$Spatial_CT_broad <- recode(
  RM_ST$Spatial_CT_broad,
  "Terminal Keratinocyte" = "Terminal keratinocytes",
  "B cells"        = "B cells",
  "Salivary"              = "Salivary gland cells"
)
RM_ST$Spatial_CT_broad <- factor(RM_ST$Spatial_CT_broad, levels = final_levels)
table(RM_ST$Spatial_CT_broad,useNA = "ifany")

broad_cols <- c(
  "Malignant cells"              = "#D41114",  # bright clinical red
  "Fibroblasts"                   = "#984EA3",  # purple (stromal)
  "Macrophages"                   = "#FF7F00",  # orange (myeloid)
  "B cells"       = "#377EB8",  # blue
  "T cells"                = "#4DAF4A",  # green
  "Dendritic cells"        = "#A65628",  # brown
  "Endothelial cells"      = "#4DD0E1",  # light cyan (distinct from all)
  "Terminal keratinocytes" = "#F781BF",  # pink (epithelial)
  "Muscle cells"           = "#999999",  # grey
  "Salivary gland cells"   = "#1B9E77"   # teal
)

p_broad_png <- umap_png(
  group_by   = "Spatial_CT_broad",
  colors     = broad_cols,
  title_text = "Broad Cell Types"
)
print(p_broad_png)

ggsave(filename = file.path(out_dir, "UMAP_BroadCellTypes.png"),
  p_broad_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)

p3 <- SpatialDimPlot(
  RM_ST,
  group.by = "Spatial_CT_broad",
  images   = "spatial.6",
  cols     = broad_cols,
  pt.size.factor = 5,
  combine  = FALSE
)[[1]]

p3 <- p3 + theme(legend.position = "none")

ggsave(filename = file.path(out_dir, "Spatial_spatial6_BroadCellTypes_noLegend.png"),
  p3,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)


unique(RM_ST$ID)
table(RM_ST$ID)

# IDs to preserve
keep_post <- c("LIB15_Post", "LIB17_Post")

# Start from original IDs
old_ids <- RM_ST$ID  

# Default: use the existing ID
new_ids <- old_ids  

# Step 1 — change ALL _Post → _OnT
new_ids <- gsub("_Post$", "_OnT", new_ids)

# Step 2 — restore the exceptions (keep them as Post)
new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]

# Assign back to object
RM_ST$Sample_ID_renamed <- new_ids


table(RM_ST$ID, RM_ST$Sample_ID_renamed,useNA = "ifany")

table(RM_ST$Sample_ID_renamed,useNA = "ifany")

RM_ST$Sample_ID_renamed <- factor(
  RM_ST$Sample_ID_renamed,
  levels = sort(unique(RM_ST$Sample_ID_renamed))
)

sample_ids <- levels(RM_ST$Sample_ID_renamed)

sample_cols <- setNames(
  hcl.colors(length(sample_ids), palette = "Dynamic"),  # or "Dark3", "Set3", "VikO"
  sample_ids
)

p_sample_png <- umap_png(
  group_by   = "Sample_ID_renamed",
  colors     = sample_cols,
  title_text = "Sample ID"
)

print(p_sample_png)



ggsave(filename = file.path(out_dir, "UMAP_SampleID.png"),
  p_sample_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)





#Figure 1def
# Where you want to save figures
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/")  # removed (publication): use out_dir + file.path()
set.seed(1234)  # for reproducibility if you ever rerun UMAP

# Keep clustering neighbors/graph as-is (dims 1:15), just make a "visual UMAP"
RM_ST <- RunUMAP(
  RM_ST,
  reduction     = "integrated.dr",
  dims          = 1:5,                  # <- fewer PCs for viz
  reduction.name = "umap.integrated.vis"
)
ElbowPlot(RM_ST)
p_vis <- DimPlot(
  RM_ST,
  reduction = "umap.integrated.vis",
  group.by  = "Pathologist_Annotation",
  raster    = F,
  pt.size   = 1
)

###############################################
## 1. Quick sanity checks
###############################################

levels(factor(RM_ST$Spatial_CT))
table(RM_ST$Spatial_CT, useNA = "ifany")

levels(factor(RM_ST$Pathologist_Annotation))
table(RM_ST$Pathologist_Annotation, useNA = "ifany")

levels(factor(RM_ST$ID))
table(RM_ST$ID, useNA = "ifany")

###############################################
## 2. Clean + order pathologist annotations
###############################################

# Make sure levels are consistent and ordered for plotting
RM_ST$Pathologist_Annotation <- factor(
  RM_ST$Pathologist_Annotation,
  levels = c("SCC","Non-SCC", "Out-of-slide")
)

# Color palette (publication-friendly)
ann_cols <- c(
  "Non-SCC"      = "#08519C",   # green
  "SCC"          = "#C51B7D",   # magenta-red (strong tumor color)
  "Out-of-slide" = "#999999"    # dark neutral gray
)


mal_cols <- c(
  "Malignant spots"     = "#FF0000",  # deep magenta (tumor)
  "Non-Malignant spots" = "darkgrey"   # green (stromal/immune)
)


#We need some help with pin-pointing the 9 OnT samples and 2 Post-Treatment Samles
#LIB15_Post => LIB15_OnT
#LIB17_Post => LIB17_OnT

DimPlot(RM_ST,reduction = "umap.integrated", group.by = "Pathologist_Annotation")
DimPlot(RM_ST,reduction = "umap.integrated", group.by = "Pathologist_Annotation")




###############################################
## 1. Define malignant Spatial_CT categories
###############################################

malignant_clusters <- c(
  "TC_Epithelial_Differentiated",
  "Neuroendocrine/Stem-like Niche",
  "LE_ECM_Remodeling",
  "LE_Pro_Inflammatory",
  "Hypoxic/Angiogenic Niche"
)

RM_ST$Consensus_Malignancy <- ifelse(
  RM_ST$Spatial_CT %in% malignant_clusters,
  "Malignant spots",
  "Non-Malignant spots"
)

# Turn into an ordered factor
RM_ST$Consensus_Malignancy <- factor(
  RM_ST$Consensus_Malignancy,
  levels = c("Malignant spots", "Non-Malignant spots")
)

table(RM_ST$Consensus_Malignancy)


############################################################
## 0. UMAP limits (run once after UMAP is computed)
############################################################

umap_df <- Embeddings(RM_ST, "umap.integrated.vis")

xlim_fixed <- range(umap_df[, 1])
ylim_fixed <- range(umap_df[, 2])

xlim_fixed
ylim_fixed

############################################################
## 1. Color palettes
############################################################

# Pathologist annotation colors
ann_cols <- c(
  "Non-SCC"      = "darkblue",   # or e.g. "#08519C" if you want hex
  "SCC"          = "#C51B7D",    # magenta-red tumor
  "Out-of-slide" = "#999999"     # grey
)

# Consensus malignancy colors with brighter red
mal_cols <- c(
  "Malignant spots"     = "#D41114",   # bright strong red
  "Non-Malignant spots" = "#999999"    # green
)

############################################################
## 2. Theme for PNG figures (no box, legend ON)
############################################################

theme_umap <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title       = element_blank(),   # remove axis titles
      axis.text        = element_blank(),   # remove axis tick labels
      axis.ticks       = element_blank(),   # remove tick marks
      axis.line        = element_blank(),   # remove axis lines
      
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.position  = "right",
      legend.key.size  = unit(0.35, "cm"),
      
      plot.title       = element_text(size = 11, face = "bold", 
                                      hjust = 0.5)
    )
}


############################################################
## 3. Standardized UMAP function for PNGs (with legend)
############################################################

umap_png <- function(group_by, colors, title_text) {
  
  DimPlot(
    RM_ST,
    reduction = "umap.integrated.vis",
    group.by  = group_by,
    pt.size   = 0.8,     # 1.5 is very chunky; 0.4–0.8 is usually ideal
    raster    = F
  ) +
    scale_color_manual(values = colors, drop = FALSE) +
    
    # Fixed coordinate range → identical panel framing
    coord_cartesian(
      xlim   = xlim_fixed,
      ylim   = ylim_fixed,
      expand = FALSE
    ) +
    
    theme_umap() +
    
    labs(
      x     = "UMAP 1",
      y     = "UMAP 2",
      color = title_text,
      title = title_text
    )
}

############################################################
## 4. Build the UMAPs
############################################################

p_path_png <- umap_png(
  group_by   = "Pathologist_Annotation",
  colors     = ann_cols,
  title_text = "Pathologist Annotation"
)
print(p_path_png)

p_consensus_png <- umap_png(
  group_by   = "Consensus_Malignancy",
  colors     = mal_cols,
  title_text = "Consensus-Based Malignancy"
)
print(p_consensus_png)

############################################################
## 5. High-resolution PNG export (identical size)
############################################################

ggsave(filename = file.path(out_dir, "UMAP_PathologistAnnotation.png"),
  p_path_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "UMAP_ConsensusMalignancy.png"),
  p_consensus_png,
  width  = 7.5,
  height = 5,
  units  = "in",
  dpi    = 900
)

p1 <- SpatialDimPlot(
  RM_ST,
  group.by = "Consensus_Malignancy",
  images = "spatial.6",
  cols = mal_cols,
  pt.size.factor = 5,
  combine = FALSE
)[[1]]

p2 <- SpatialDimPlot(
  RM_ST,
  group.by = "Pathologist_Annotation",
  images = "spatial.6",
  cols = ann_cols,
  pt.size.factor = 5,
  combine = FALSE
)[[1]]


p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")


ggsave(filename = file.path(out_dir, "Spatial_spatial6_ConsensusMalignancy_noLegend.png"),
  p1,
  width = 7,
  height = 5,
  units = "in",
  dpi = 900
)

ggsave(filename = file.path(out_dir, "Spatial_spatial6_PathologistAnnotation_noLegend.png"),
  p2,
  width = 7,
  height = 5,
  units = "in",
  dpi = 900
)


#Figures
malignant_clusters <- c(
  "TC_Epithelial_Differentiated",
  "Neuroendocrine/Stem-like Niche",
  "LE_ECM_Remodeling",
  "LE_Pro_Inflammatory",
  "Hypoxic/Angiogenic Niche"
)

caf_clusters <- c(
  "apCAFs",
  "mCAFs",
  "myCAFs",
  "Hypoxic/Invasion-associated CAFs (high activation)",
  "Hypoxic/Invasion-associated CAFs (transitional)",
  "Erythroid/Platelet-interacting CAFs"
)

tam_clusters <- c(
  "M1-like/Inflammatory TAMs",
  "M2-like/ECM-remodeling TAMs",
  "tTAMs"
)

bplasma_clusters <- c(
  "Mature B cell/Plasma cell",
  "Mucosal Plasma B cells",
  "tBcells"
)


RM_ST$Spatial_CT_broad <- case_when(
  RM_ST$Spatial_CT %in% malignant_clusters ~ "Malignant cells",
  RM_ST$Spatial_CT %in% caf_clusters       ~ "Fibroblasts",
  RM_ST$Spatial_CT %in% tam_clusters       ~ "Macrophages",
  RM_ST$Spatial_CT %in% bplasma_clusters   ~ "B cells",
  
  # everything else retains original identity
  TRUE ~ RM_ST$Spatial_CT
)

final_levels <- c(
  "Malignant cells",
  "Fibroblasts",
  "Macrophages",
  "B cells",
  "T cells",
  "Dendritic cells",
  "Endothelial cells",
  "Terminal keratinocytes",
  "Muscle cells",
  "Salivary gland cells"
)

RM_ST$Spatial_CT_broad <- recode(
  RM_ST$Spatial_CT_broad,
  "Terminal Keratinocyte" = "Terminal keratinocytes",
  "B cells"        = "B cells",
  "Salivary"              = "Salivary gland cells"
)
RM_ST$Spatial_CT_broad <- factor(RM_ST$Spatial_CT_broad, levels = final_levels)
table(RM_ST$Spatial_CT_broad,useNA = "ifany")

broad_cols <- c(
  "Malignant cells"              = "#D41114",  # bright clinical red
  "Fibroblasts"                   = "#984EA3",  # purple (stromal)
  "Macrophages"                   = "#FF7F00",  # orange (myeloid)
  "B cells"       = "#377EB8",  # blue
  "T cells"                = "#4DAF4A",  # green
  "Dendritic cells"        = "#A65628",  # brown
  "Endothelial cells"      = "#E41A1C",  # red-crimson (distinct from malignant)
  "Terminal keratinocytes" = "#F781BF",  # pink (epithelial)
  "Muscle cells"           = "#999999",  # grey
  "Salivary gland cells"   = "#1B9E77"   # teal
)

p_broad_png <- umap_png(
  group_by   = "Spatial_CT_broad",
  colors     = broad_cols,
  title_text = "Broad Cell Types"
)
print(p_broad_png)

ggsave(filename = file.path(out_dir, "UMAP_BroadCellTypes.png"),
  p_broad_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)

p3 <- SpatialDimPlot(
  RM_ST,
  group.by = "Spatial_CT_broad",
  images   = "spatial.6",
  cols     = broad_cols,
  pt.size.factor = 5,
  combine  = FALSE
)[[1]]

p3 <- p3 + theme(legend.position = "none")

ggsave(filename = file.path(out_dir, "Spatial_spatial6_BroadCellTypes_noLegend.png"),
  p3,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)


unique(RM_ST$ID)
table(RM_ST$ID)

# IDs to preserve
keep_post <- c("LIB15_Post", "LIB17_Post")

# Start from original IDs
old_ids <- RM_ST$ID  

# Default: use the existing ID
new_ids <- old_ids  

# Step 1 — change ALL _Post → _OnT
new_ids <- gsub("_Post$", "_OnT", new_ids)

# Step 2 — restore the exceptions (keep them as Post)
new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]

# Assign back to object
RM_ST$Sample_ID_renamed <- new_ids


table(RM_ST$ID, RM_ST$Sample_ID_renamed,useNA = "ifany")

table(RM_ST$Sample_ID_renamed,useNA = "ifany")

RM_ST$Sample_ID_renamed <- factor(
  RM_ST$Sample_ID_renamed,
  levels = sort(unique(RM_ST$Sample_ID_renamed))
)

sample_ids <- levels(RM_ST$Sample_ID_renamed)

sample_cols <- setNames(
  hcl.colors(length(sample_ids), palette = "Dynamic"),  # or "Dark3", "Set3", "VikO"
  sample_ids
)

p_sample_png <- umap_png(
  group_by   = "Sample_ID_renamed",
  colors     = sample_cols,
  title_text = "Sample ID"
)

print(p_sample_png)



ggsave(filename = file.path(out_dir, "UMAP_SampleID.png"),
  p_sample_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)


###############################################
## EXTENSION: UMAPs by Patient and Timepoint
## - Uses your existing umap_png(), theme_umap(),
##   xlim_fixed/ylim_fixed, and umap.integrated.vis
###############################################



suppressPackageStartupMessages({
  library(stringr)
})

# 0) Prefer Sample_ID_renamed else ID
id_col <- if ("Sample_ID_renamed" %in% colnames(RM_ST@meta.data)) "Sample_ID_renamed" else "ID"

# 1) Pull as a NAMED atomic vector (names = cell barcodes)
id_vec <- FetchData(RM_ST, vars = id_col)[, 1]
id_vec <- as.character(id_vec)
names(id_vec) <- colnames(RM_ST)

# Safety check (should be TRUE)
stopifnot(identical(names(id_vec), colnames(RM_ST)))

# 2) Parse Patient_ID + Timepoint
patient_vec <- str_replace(id_vec, "_[^_]+$", "")
timepoint_vec <- str_extract(id_vec, "[^_]+$")

# 3) Optional normalize timepoint labels
timepoint_vec <- recode(
  timepoint_vec,
  "BL"   = "Pre",
  "Pre"  = "Pre",
  "OnT"  = "On-treatment",
  "Post" = "Post-treatment",
  .default = timepoint_vec
)

# 4) Add back to Seurat meta.data with guaranteed alignment
RM_ST[["Patient_ID"]] <- patient_vec
RM_ST[["Timepoint"]]  <- timepoint_vec

# 5) Make factors (stable ordering)
RM_ST$Patient_ID <- factor(RM_ST$Patient_ID, levels = sort(unique(RM_ST$Patient_ID)))
RM_ST$Timepoint  <- factor(RM_ST$Timepoint, levels = c("Pre","On-treatment","Post-treatment"))
RM_ST$Timepoint  <- droplevels(RM_ST$Timepoint)

# Sanity
table(RM_ST$Patient_ID, useNA = "ifany")
table(RM_ST$Timepoint, useNA = "ifany")

# ----------------------------
# 2) Palettes
#    Patient palette: dynamic, many colors
#    Timepoint palette: fixed, publication-friendly
# ----------------------------
patient_levels <- levels(RM_ST$Patient_ID)

patient_cols <- setNames(
  hcl.colors(length(patient_levels), palette = "Dynamic"),
  patient_levels
)

timepoint_levels <- levels(RM_ST$Timepoint)

# Feel free to tweak these 3 colors to match your paper style
timepoint_cols <- c(
  "Pre"        = "#1F78B4",  # blue
  "On-treatment"    = "#FF7F00",  # orange
  "Post-treatment"  = "#33A02C"   # green
)
# Keep only what exists
timepoint_cols <- timepoint_cols[timepoint_levels]

# ----------------------------
# 3) Build the 2 new UMAPs using your wrapper
#    (expects umap_png() already defined above)
# ----------------------------
p_patient_png <- umap_png(
  group_by   = "Patient_ID",
  colors     = patient_cols,
  title_text = "Patient"
)
print(p_patient_png)

p_timepoint_png <- umap_png(
  group_by   = "Timepoint",
  colors     = timepoint_cols,
  title_text = "Timepoint"
)
print(p_timepoint_png)

# ----------------------------
# 4) Export high-res PNGs (same size/dpi as your others)
# ----------------------------
ggsave(filename = file.path(out_dir, "UMAP_Patient.png"),
  p_patient_png,
  width  = 6.5,
  height = 5,
  units  = "in",
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "UMAP_Timepoint.png"),
  p_timepoint_png,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)

levels(factor(RM_ST$Timepoint))






