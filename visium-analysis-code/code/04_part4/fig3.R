# ============================================================
# Figure 3 — Non-malignant compartment panels
# Publication polished script (paths parameterized, no setwd)
# ============================================================


# ---- shared helpers ----
# Assumes this script lives in the same folder as 00_common.R.
# If not, change the path below (e.g., source("code/04_part4/00_common.R")).
source("00_common.R")
theme_set(theme_pub())

suppressPackageStartupMessages({
  library(Nebulosa)
  library(SCpubr)
  library(Seurat)
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
  library(phangorn)
  library(treeio)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------
input_rds_1 <- "part2/Fibroblast_v2.RDS"
input_rds_2 <- "part2/Macrophage.RDS"
input_rds_3 <- "part2/B_cells.RDS"

# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "fig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


###############################################
## 0. Subset Non-Malignant Spots
###############################################

# 1) Make a copy without images
RM_ST_noimg <- RM_ST
RM_ST_noimg@images <- list()   # drop Visium images so subset doesn't validate them

# Keep all broad cell types except "Malignant"
RM_ST_nonmalig <- subset(
  RM_ST_noimg,
  subset = Spatial_CT_broad != "Malignant cells"
)

table(RM_ST_nonmalig$Spatial_CT_broad, useNA = "ifany")
table(RM_ST_noimg$Spatial_CT_broad, useNA = "ifany")

# Drop unused factor levels for cleaner plotting
RM_ST_nonmalig$Spatial_CT_broad        <- droplevels(RM_ST_nonmalig$Spatial_CT_broad)
RM_ST_nonmalig$Sample_ID_renamed       <- droplevels(RM_ST_nonmalig$Sample_ID_renamed)
RM_ST_nonmalig$Pathologist_Annotation  <- droplevels(RM_ST_nonmalig$Pathologist_Annotation)

table(RM_ST_nonmalig$Spatial_CT_broad, useNA = "ifany")

###############################################
## 1. UMAP for Non-Malignant Spots
###############################################

set.seed(1234)

RM_ST_nonmalig <- RunUMAP(
  RM_ST_nonmalig,
  reduction      = "integrated.dr",   # same PCA space as full object
  dims           = 1:5,               # consistent with Fig 1
  reduction.name = "umap.integrated.nonmalig"
)

# Fix UMAP limits for consistent framing across panels
umap_nonmalig_df <- Embeddings(RM_ST_nonmalig, "umap.integrated.nonmalig")
xlim_nonmalig    <- range(umap_nonmalig_df[, 1])
ylim_nonmalig    <- range(umap_nonmalig_df[, 2])

xlim_nonmalig
ylim_nonmalig


###############################################
## 2. Colors specific to Non-Malignant subset
###############################################

# Broad cell types, excluding malignant
broad_cols_nonmalig <- broad_cols[
  names(broad_cols) != "Malignant cells"
]

# Sample colors for only the samples that remain after subsetting
sample_ids_nonmalig <- levels(RM_ST_nonmalig$Sample_ID_renamed)
sample_cols_nonmalig <- sample_cols[sample_ids_nonmalig]

# Pathologist colors can be reused directly;
# levels may have dropped (e.g. maybe no "SCC" in this subset)
ann_cols_nonmalig <- ann_cols[names(ann_cols) %in% levels(RM_ST_nonmalig$Pathologist_Annotation)]
ann_cols_nonmalig



###############################################
## 3. Standardized UMAP function (Non-Malignant)
###############################################

umap_nonmalig_png <- function(group_by, colors, title_text) {
  DimPlot(
    RM_ST_nonmalig,
    reduction = "umap.integrated.nonmalig",
    group.by  = group_by,
    pt.size   = 0.8,
    raster    = FALSE
  ) +
    scale_color_manual(values = colors, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_nonmalig,
      ylim   = ylim_nonmalig,
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

###############################################
## 4. Generate Non-Malignant UMAP Panels
###############################################

# (a) Broad non-malignant compartments
p_nm_broad <- umap_nonmalig_png(
  group_by   = "Spatial_CT_broad",
  colors     = broad_cols_nonmalig,
  title_text = "Broad Cell Types (Non-Malignant)"
)
print(p_nm_broad)

# (b) Sample ID
p_nm_sample <- umap_nonmalig_png(
  group_by   = "Sample_ID_renamed",
  colors     = sample_cols_nonmalig,
  title_text = "Sample ID (Non-Malignant Spots)"
)
print(p_nm_sample)


###############################################
## 5. Export Non-Malignant UMAPs
###############################################
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "UMAP_NonMalignant_BroadCellTypes.png"),
  p_nm_broad,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "UMAP_NonMalignant_SampleID.png"),
  p_nm_sample,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)


###############################################
## 6. Marker groups for Non-Malignant density plots
###############################################

# Group 1 – Fibroblast / CAF markers
fibro_genes <- c("FAP", "COL1A1", "PDGFRB")

# Group 2 – B / Plasma markers (narrowed to 3)
bplas_genes <- c("CD19", "MS4A1", "CD38")

# Group 3 – Macrophage / Myeloid markers
mac_genes <- c("CD68", "CD14", "CSF1R")

# Group 4 – T cell markers
tcell_genes <- c("CD3D", "CD4", "CD8A")



DefaultAssay(RM_ST_nonmalig) <- "Spatial"   # or "SCT"/"RNA" depending where these live

gene_list_all <- c(fibro_genes, bplas_genes, mac_genes, tcell_genes)
gene_list_all[!gene_list_all %in% rownames(GetAssayData(RM_ST_nonmalig))]
# ideally returns character(0)


###############################################
## 7. Nebulosa density plots on RM_ST_nonmalig
###############################################


# make sure the assay is set (already done above if you ran that line)
# DefaultAssay(RM_ST_nonmalig) <- "Spatial"

nm_neb_plot_single_titles <- function(features) {
  SCpubr::do_NebulosaPlot(
    sample               = RM_ST_nonmalig,
    features             = features,
    reduction            = "umap.integrated.nonmalig",
    use_viridis          = FALSE,
    sequential.palette   = "YlGnBu",
    sequential.direction = 1,
    legend.position      = "bottom",
    plot.axes            = FALSE,
    font.size            = 25,
  ) 
}



# Build panels
p_nm_fibro <- nm_neb_plot_single_titles(fibro_genes)
p_nm_bplas  <- nm_neb_plot_single_titles(bplas_genes)
p_nm_mac    <- nm_neb_plot_single_titles(mac_genes)
p_nm_tcells <- nm_neb_plot_single_titles(tcell_genes)



p_nm_fibro <- nm_neb_plot_single_titles(fibro_genes) +
  plot_annotation(title = "Fibroblast Markers") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_bplas <- nm_neb_plot_single_titles(bplas_genes) +
  plot_annotation(title = "B cell Markers") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_mac <- nm_neb_plot_single_titles(mac_genes) +
  plot_annotation(title = "Macrophage Markers") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_tcells <- nm_neb_plot_single_titles(tcell_genes) +
  plot_annotation(title = "T cell Markers") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "NM_Fibroblast_CAF_markers.png"),
  p_nm_fibro,
  width  = 13,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_B_Plasma_markers.png"),
  p_nm_bplas,
  width  = 13,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_Macrophage_markers.png"),
  p_nm_mac,
  width  = 13,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_Tcell_markers.png"),
  p_nm_tcells,
  width  = 13,
  height = 6,
  dpi    = 900
)




#Fibroblast

Fibroblast <- readRDS(input_rds_1)

DimPlot(Fibroblast, reduction = "umap.cca_Fibroblasts_v2", pt.size = 1.5, group.by = "Fibroblast_SubCluster")

table(Fibroblast$Fibroblast_SubCluster, useNA = "ifany")

fibro_remap <- c(
  "Hypoxic/Invasion-associated CAFs (transitional)" = "cluster 2",
  "mCAFs"                                           = "cluster 6",
  "apCAFs"                                          = "cluster 3",
  "Hypoxic/Invasion-associated CAFs (high activation)" = "cluster 5",
  "myCAFs"                                          = "cluster 1",
  "Erythroid/Platelet-interacting CAFs"             = "cluster 4"
)

# This creates a *named* vector; we drop the names with unname()
Fibroblast@meta.data$Fibroblast_Cluster <- unname(
  fibro_remap[ Fibroblast$Fibroblast_SubCluster ]
)

Fibroblast@meta.data$Fibroblast_Cluster <- factor(
  Fibroblast@meta.data$Fibroblast_Cluster,
  levels = paste0("cluster ", 1:6)
)

table(Fibroblast@meta.data$Fibroblast_Cluster, useNA = "ifany")




# Check reduction name
Reductions(Fibroblast)
# assuming: "umap.cca_Fibro" already exists

umap_fibro <- "umap.cca_Fibroblasts_v2"

# Get UMAP coordinates and padded limits
umap_fibro_df <- Embeddings(Fibroblast, umap_fibro)
xlim_fibro <- range(umap_fibro_df[, 1])
ylim_fibro <- range(umap_fibro_df[, 2])

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_fibro_pad <- pad_limits(xlim_fibro, 0.03)
ylim_fibro_pad <- pad_limits(ylim_fibro, 0.03)

# Cluster levels & colors (Okabe–Ito style)
cluster_levels_fibro <- levels(Fibroblast$Fibroblast_Cluster)

fibro_cluster_colors <- setNames(
  c( "#0072B2", "#009E73", "#CC79A7", "#F0E442", "#D55E00" ,"#56B4E9"),
  cluster_levels_fibro
)

# UMAP theme (same style as malignant)
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
  group.by  = "Fibroblast_Cluster",
  pt.size   = 0.8,
  raster    = FALSE
) +
  scale_color_manual(values = fibro_cluster_colors, drop = FALSE) +
  coord_cartesian(
    xlim   = xlim_fibro_pad,
    ylim   = ylim_fibro_pad,
    expand = FALSE
  ) +
  theme_umap_fibro() +
  labs(
    x     = NULL,
    y     = NULL,
    color = "Fibroblast cluster",
    title = "Fibroblast UMAP (cluster 1–6)"
  )

fibro_umap
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Fig3cd/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "UMAP_Fibroblast_Clusters.png"),
  fibro_umap,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)




# Use SCT or RNA depending on what you used originally
DefaultAssay(Fibroblast) <- "SCT"

# Set identities to the reverted labels
Idents(Fibroblast) <- "Fibroblast_Cluster"

# Build cluster tree using PCs (same dims as you used for Cancer_v3)
Fibroblast <- BuildClusterTree(
  object  = Fibroblast,
  dims    = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

fibro_tree <- Fibroblast@tools$BuildClusterTree  # hclust
fibro_phy  <- as.phylo(fibro_tree)

fibro_cluster_levels <- levels(Fibroblast$Fibroblast_Cluster)

# Convert tree to tibble and annotate internal nodes by majority tip cluster
fibro_tree_dat <- as_tibble(fibro_phy)

get_leaf_descendants <- function(node, phy) {
  phangorn::Descendants(phy, node, type = "tips")[[1]]
}

fibro_tree_dat$cluster <- sapply(fibro_tree_dat$node, function(node) {
  tip_ids <- get_leaf_descendants(node, fibro_phy)
  if (length(tip_ids) == 0) return(NA_character_)
  
  tip_labels <- fibro_phy$tip.label[tip_ids]
  tip_labels <- tip_labels[tip_labels %in% fibro_cluster_levels]
  if (length(tip_labels) == 0) return(NA_character_)
  
  tbl <- table(tip_labels)
  if (length(tbl) == 1) {
    return(names(tbl)[1])
  } else {
    return(names(tbl)[which.max(tbl)])  # majority cluster
  }
})

fibro_tree_dat$cluster <- factor(
  fibro_tree_dat$cluster,
  levels = fibro_cluster_levels
)

# reuse fibro_cluster_colors defined above
fibro_tree_dat$pretty_label <- fibro_tree_dat$label

p_fibro_tree <- ggtree(fibro_phy, aes(color = cluster)) %<+% fibro_tree_dat +
  geom_tree(size = 1.5) +
  ggtree::geom_point2(aes(subset = !is.na(cluster)), size = 5) +
  geom_tiplab(
    aes(label = pretty_label, color = cluster),
    offset   = 0.5,
    size     = 6,
    fontface = "bold"
  ) +
  scale_color_manual(
    values   = fibro_cluster_colors,
    na.value = "grey80",
    name     = "Fibroblast cluster"
  ) +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  theme_tree() +
  theme(
    plot.margin     = unit(c(0.2, 2.5, 0.2, 0.2), "cm"),
    legend.position = "none"
  )

p_fibro_tree





ggsave(filename = file.path(out_dir, "Tree_Fibroblast_Clusters_clean.png"),
  p_fibro_tree,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)




#Macrophage

Macrophage <- readRDS(input_rds_2)

DimPlot(Macrophage, reduction = "umap.cca_macro", pt.size = 1.5, group.by = "Macrophage_Identity")

table(Macrophage$Macrophage_Identity, useNA = "ifany")

macro_remap <- c(
  "M2-like/ECM-remodeling TAMs" = "cluster 1",
  "M1-like/Inflammatory TAMs"   = "cluster 2",
  "tTAMs"                       = "cluster 3"
)


# This creates a *named* vector; we drop the names with unname()
Macrophage@meta.data$Macrophage_Cluster <- unname(
  macro_remap[ Macrophage$Macrophage_Identity ]
)

Macrophage@meta.data$Macrophage_Cluster <- factor(
  Macrophage@meta.data$Macrophage_Cluster,
  levels = paste0("cluster ", 1:3)
)

table(Macrophage@meta.data$Macrophage_Cluster, useNA = "ifany")




# Check reduction name
Reductions(Macrophage)
# assuming: "umap.cca_Fibro" already exists

umap_macro <- "umap.cca_macro"

# Get UMAP coordinates and padded limits
umap_macro_df <- Embeddings(Macrophage, umap_macro)
xlim_macro <- range(umap_macro_df[, 1])
ylim_macro <- range(umap_macro_df[, 2])

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_macro_pad <- pad_limits(xlim_macro, 0.03)
ylim_macro_pad <- pad_limits(ylim_macro, 0.03)

# Cluster levels & colors 
cluster_levels_macro <- levels(Macrophage$Macrophage_Cluster)

## 1. Macrophage-specific colors (no overlap with fibro Okabe–Ito)
macro_cluster_colors <- c(
  "cluster 1" = "#CC6677",  # reddish rose
  "cluster 2" = "#882255",  # deep wine
  "cluster 3" = "#44AA99"   # teal
)
macro_cluster_colors

# UMAP theme (matching your previous style)
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
  group.by  = "Macrophage_Cluster",
  pt.size   = 0.8,
  raster    = FALSE
) +
  scale_color_manual(
    values = macro_cluster_colors,
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
    color = "Macrophage cluster",
    title = "Macrophage UMAP (cluster 1-3)"
  )

macro_umap
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Fig3cd/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "UMAP_Macrophage_Clusters.png"),
  macro_umap,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)



## Use same assay as in fibro / cancer trees
DefaultAssay(Macrophage) <- "SCT"   # or "RNA" if that’s what you used

Idents(Macrophage) <- "Macrophage_Cluster"

Macrophage <- BuildClusterTree(
  object  = Macrophage,
  dims    = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

macro_tree <- Macrophage@tools$BuildClusterTree  # hclust
macro_phy  <- as.phylo(macro_tree)

macro_cluster_levels <- levels(Macrophage$Macrophage_Cluster)

macro_tree_dat <- as_tibble(macro_phy)

get_leaf_descendants <- function(node, phy) {
  phangorn::Descendants(phy, node, type = "tips")[[1]]
}

macro_tree_dat$cluster <- sapply(macro_tree_dat$node, function(node) {
  tip_ids <- get_leaf_descendants(node, macro_phy)
  if (length(tip_ids) == 0) return(NA_character_)
  
  tip_labels <- macro_phy$tip.label[tip_ids]
  tip_labels <- tip_labels[tip_labels %in% macro_cluster_levels]
  if (length(tip_labels) == 0) return(NA_character_)
  
  tbl <- table(tip_labels)
  if (length(tbl) == 1) {
    return(names(tbl)[1])
  } else {
    return(names(tbl)[which.max(tbl)])  # majority cluster
  }
})

macro_tree_dat$cluster <- factor(
  macro_tree_dat$cluster,
  levels = macro_cluster_levels
)

macro_tree_dat$pretty_label <- macro_tree_dat$label

p_macro_tree <- ggtree(macro_phy, aes(color = cluster)) %<+% macro_tree_dat +
  geom_tree(size = 1.5) +
  ggtree::geom_point2(aes(subset = !is.na(cluster)), size = 5) +
  geom_tiplab(
    aes(label = pretty_label, color = cluster),
    offset   = 0.5,
    size     = 6,
    fontface = "bold"
  ) +
  scale_color_manual(
    values   = macro_cluster_colors,
    na.value = "grey80",
    name     = "Macrophage cluster"
  ) +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  theme_tree() +
  theme(
    plot.margin     = unit(c(0.2, 2.5, 0.2, 0.2), "cm"),
    legend.position = "none"
  )

p_macro_tree



ggsave(filename = file.path(out_dir, "Tree_Macrophage_Clusters_clean.png"),
  p_macro_tree,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)



###############################################
###############################################

Cluster0 <- c("CXCL9", "MMP9")
Cluster1 <- c("IL1B", "CXCL2")
Cluster2 <- c("S100A8", "S100A9")

###############################################
## 7. Nebulosa density plots on RM_ST_nonmalig
###############################################


# make sure the assay is set (already done above if you ran that line)
DefaultAssay(Macrophage) <- "Spatial"

nm_neb_plot_single_titles <- function(features) {
  SCpubr::do_NebulosaPlot(
    sample               = Macrophage,
    features             = features,
    reduction            = "umap.cca_macro",
    use_viridis          = FALSE,
    sequential.palette   = "YlGnBu",
    sequential.direction = 1,
    legend.position      = "bottom",
    plot.axes            = FALSE,
    font.size            = 25,
  ) 
}


nm_neb_plot_single_titles <- function(features) {
  SCpubr::do_NebulosaPlot(
    sample               = Fibroblast,
    features             = features,
    reduction            = "umap.cca_Fibroblasts_v2",
    use_viridis          = FALSE,
    sequential.palette   = "YlGnBu",
    sequential.direction = 1,
    legend.position      = "bottom",
    plot.axes            = FALSE,
    font.size            = 25,
  ) 
}

Fibro_gene_cluster1 <- c("LRRC15", "MMP11")
Fibro_gene_cluster2 <- c("CCL19", "CXCR4")
Fibro_gene_cluster3 <- c("IL1B", "CSF3")
Fibro_gene_cluster4 <- c("MYH11", "ACTG2")

p_nm_fibro_C1 <- nm_neb_plot_single_titles(Fibro_gene_cluster1)
p_nm_fibro_C1 <- nm_neb_plot_single_titles(Fibro_gene_cluster1) +
  plot_annotation(title = "ecm-myCAF") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_fibro_C2 <- nm_neb_plot_single_titles(Fibro_gene_cluster2)
p_nm_fibro_C2 <- nm_neb_plot_single_titles(Fibro_gene_cluster2) +
  plot_annotation(title = "IFNg-iCAF") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_fibro_C3 <- nm_neb_plot_single_titles(Fibro_gene_cluster3)
p_nm_fibro_C3 <- nm_neb_plot_single_titles(Fibro_gene_cluster3) +
  plot_annotation(title = "IL-iCAF") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_fibro_C4 <- nm_neb_plot_single_titles(Fibro_gene_cluster4)
p_nm_fibro_C4 <- nm_neb_plot_single_titles(Fibro_gene_cluster4) +
  plot_annotation(title = "acto-myCAF") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Supp_Fig5/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "NM_fibroblast_Cluster1.png"),
  p_nm_fibro_C1,
  width  = 8.5,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_fibroblast_Cluster2.png"),
  p_nm_fibro_C2,
  width  = 8.5,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_fibroblast_Cluster3.png"),
  p_nm_fibro_C3,
  width  = 8.5,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_fibroblast_Cluster4.png"),
  p_nm_fibro_C4,
  width  = 8.5,
  height = 6,
  dpi    = 900
)



###############################################

# Build panels
p_nm_Cluster0 <- nm_neb_plot_single_titles(Cluster0)
p_nm_Cluster1  <- nm_neb_plot_single_titles(Cluster1)
p_nm_Cluster2  <- nm_neb_plot_single_titles(Cluster2)

p_nm_Cluster0 <- nm_neb_plot_single_titles(Cluster0) +
  plot_annotation(title = "IFN-TAM") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_Cluster1 <- nm_neb_plot_single_titles(Cluster1) +
  plot_annotation(title = "pro-inflammatory TAM") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

p_nm_Cluster2 <- nm_neb_plot_single_titles(Cluster2) +
  plot_annotation(title = "monocyte") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Supp_Fig5/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "NM_Macrophage_Cluster0.png"),
  p_nm_Cluster0,
  width  = 8.5,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_Macrophage_Cluster1.png"),
  p_nm_Cluster1,
  width  = 8.5,
  height = 6,
  dpi    = 900
)

ggsave(filename = file.path(out_dir, "NM_Macrophage_Cluster2.png"),
  p_nm_Cluster2,
  width  = 8.5,
  height = 6,
  dpi    = 900
)

#B cell

B_cells <- readRDS(input_rds_3)

DimPlot(B_cells, reduction = "umap.cca_bcell", pt.size = 1.5, group.by = "Bcell_Identity")

table(B_cells$Bcell_Identity, useNA = "ifany")

# B cell remap (already defined)
bcell_remap <- c(
  "tBcells"                   = "Bcell_C0",
  "Mature B cell/Plasma cell" = "Bcell_C1",
  "Mucosal Plasma B cells"    = "Bcell_C2"
)

# Correctly map Bcell_Identity -> Bcell_C0–C2
B_cells@meta.data$Bcell_Cluster <- unname(
  bcell_remap[ B_cells$Bcell_Identity ]
)

B_cells@meta.data$Bcell_Cluster <- factor(
  B_cells@meta.data$Bcell_Cluster,
  levels = paste0("Bcell_C", 0:2)
)

table(B_cells@meta.data$Bcell_Cluster, useNA = "ifany")
# should now show counts for C0, C1, C2 and 0 NAs



# Check reduction name
Reductions(B_cells)
# assuming: "umap.cca_Fibro" already exists

umap_bcell <- "umap.cca_bcell"
umap_bcell_df <- Embeddings(B_cells, umap_bcell)
xlim_bcell <- range(umap_bcell_df[, 1])
ylim_bcell <- range(umap_bcell_df[, 2])

pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_bcell_pad <- pad_limits(xlim_bcell, 0.03)
ylim_bcell_pad <- pad_limits(ylim_bcell, 0.03)


# Cluster levels & colors 
cluster_levels_bcell <- levels(B_cells$Bcell_Cluster)

bcell_cluster_colors <- c(
  "Bcell_C0" = "#332288",  # deep indigo
  "Bcell_C1" = "#E69F00",  # ochre/orange
  "Bcell_C2" = "#999999"   # neutral grey
)

bcell_cluster_colors

# UMAP theme (matching your previous style)
theme_umap_bcell <- function() {
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

bcell_umap <- DimPlot(
  B_cells,
  reduction = umap_bcell,
  group.by  = "Bcell_Cluster",
  pt.size   = 0.8,
  raster    = FALSE
) +
  scale_color_manual(
    values = bcell_cluster_colors,
    drop   = FALSE
  ) +
  coord_cartesian(
    xlim   = xlim_bcell_pad,
    ylim   = ylim_bcell_pad,
    expand = FALSE
  ) +
  theme_umap_bcell() +
  labs(
    x     = NULL,
    y     = NULL,
    color = "B cell cluster",
    title = "B-cell UMAP (Bcell_C0–C2)"
  )

bcell_umap
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/")  # removed (publication): use out_dir + file.path()
ggsave(filename = file.path(out_dir, "UMAP_Bcell_Clusters.png"),
  bcell_umap,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)



DefaultAssay(B_cells) <- "SCT"   # or "RNA" if that’s what you used

Idents(B_cells) <- "Bcell_Cluster"

B_cells <- BuildClusterTree(
  object  = B_cells,
  dims    = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

bcell_tree <- B_cells@tools$BuildClusterTree  # hclust
bcell_phy  <- as.phylo(bcell_tree)

bcell_cluster_levels <- levels(B_cells$Bcell_Cluster)

bcell_tree_dat <- as_tibble(bcell_phy)

get_leaf_descendants <- function(node, phy) {
  phangorn::Descendants(phy, node, type = "tips")[[1]]
}

bcell_tree_dat$cluster <- sapply(bcell_tree_dat$node, function(node) {
  tip_ids <- get_leaf_descendants(node, bcell_phy)
  if (length(tip_ids) == 0) return(NA_character_)
  
  tip_labels <- bcell_phy$tip.label[tip_ids]
  tip_labels <- tip_labels[tip_labels %in% bcell_cluster_levels]
  if (length(tip_labels) == 0) return(NA_character_)
  
  tbl <- table(tip_labels)
  if (length(tbl) == 1) {
    return(names(tbl)[1])
  } else {
    return(names(tbl)[which.max(tbl)])  # majority cluster
  }
})

bcell_tree_dat$cluster <- factor(
  bcell_tree_dat$cluster,
  levels = bcell_cluster_levels
)

bcell_tree_dat$pretty_label <- bcell_tree_dat$label

p_bcell_tree <- ggtree(bcell_phy, aes(color = cluster)) %<+% bcell_tree_dat +
  geom_tree(size = 1.2) +
  ggtree::geom_point2(aes(subset = !is.na(cluster)), size = 2.5) +
  geom_tiplab(
    aes(label = pretty_label, color = cluster),
    offset   = 0.5,
    size     = 3,
    fontface = "bold"
  ) +
  scale_color_manual(
    values   = bcell_cluster_colors,
    na.value = "grey80",
    name     = "B cell cluster"
  ) +
  scale_y_reverse() +
  coord_cartesian(clip = "off") +
  theme_tree() +
  theme(
    plot.margin     = unit(c(0.2, 2.5, 0.2, 0.2), "cm"),
    legend.position = "none"
  )

p_bcell_tree


ggsave(filename = file.path(out_dir, "Tree_Bcell_Clusters_clean.png"),
  p_bcell_tree,
  width  = 7,
  height = 5,
  units  = "in",
  dpi    = 900
)