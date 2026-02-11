# ============================================================
# Figure 2 — Malignant signatures & subclustering panels
# Publication polished script (paths parameterized, no setwd)
# ============================================================


# ---- shared helpers ----
# Assumes this script lives in the same folder as 00_common.R.
# If not, change the path below (e.g., source("code/04_part4/00_common.R")).
source("00_common.R")
theme_set(theme_pub())

suppressPackageStartupMessages({
  library(AUCell)
  library(BiocGenerics)
  library(BiocManager)
  library(Matrix)
  library(Nebulosa)
  library(RColorBrewer)
  library(SCpubr)
  library(Seurat)
  library(ape)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(grid)
  library(patchwork)
  library(phangorn)
  library(pheatmap)
  library(presto)
  library(readxl)
  library(tibble)
  library(tidyr)
  library(treeio)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------
input_rds_1 <- "/part2/Cancer.rds"

# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "fig2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


#Figure 2abc (Updating Cluster #)

Cancer_v3 <- readRDS(input_rds_1)


# IDs to preserve
keep_post <- c("LIB15_Post", "LIB17_Post")

# Start from original IDs
old_ids <- Cancer_v3$ID  

# Default: use the existing ID
new_ids <- old_ids  

# Step 1 — change ALL _Post → _OnT
new_ids <- gsub("_Post$", "_OnT", new_ids)

# Step 2 — restore the exceptions (keep them as Post)
new_ids[old_ids %in% keep_post] <- old_ids[old_ids %in% keep_post]

# Assign back to object
Cancer_v3$Sample_ID_renamed <- new_ids


table(Cancer_v3$ID, Cancer_v3$Sample_ID_renamed,useNA = "ifany")

table(Cancer_v3$Sample_ID_renamed,useNA = "ifany")

Cancer_v3$Sample_ID_renamed <- factor(
  Cancer_v3$Sample_ID_renamed,
  levels = sort(unique(Cancer_v3$Sample_ID_renamed))
)

levels(factor(Cancer_v3$Cancer_Unsupervised_Clustering_v3))
levels(factor(Cancer_v3$Final_Malignant_SubCluster))
levels(factor(Cancer_v3$Sample_ID_renamed))


rename_map_cancer <- c(
  "Cancer_C0" = "cluster 5",
  "Cancer_C1" = "cluster 1",
  "Cancer_C2" = "cluster 2",
  "Cancer_C3" = "cluster 3",
  "Cancer_C4" = "cluster 4"
)


Cancer_v3$Cancer_Unsupervised_Clustering_v3_renamed <- 
  dplyr::recode(Cancer_v3$Cancer_Unsupervised_Clustering_v3, !!!rename_map_cancer)

Cancer_v3$Cancer_Unsupervised_Clustering_v3_renamed <- factor(
  Cancer_v3$Cancer_Unsupervised_Clustering_v3_renamed,
  levels = paste("cluster", 1:5)
)

table(Cancer_v3$Cancer_Unsupervised_Clustering_v3,
      Cancer_v3$Cancer_Unsupervised_Clustering_v3_renamed)


Cancer_v3$Unsupervised_Clusters <- Cancer_v3$Cancer_Unsupervised_Clustering_v3_renamed



hpv_levels <- levels(factor(Cancer_v3$hpv_status))
hpv_levels

# Start from the original column
table(Cancer_v3$hpv_status)

Cancer_v3$HPV_status_collapsed <- as.character(Cancer_v3$hpv_status)

# Treat "Non-tested" as "Negative"
Cancer_v3$HPV_status_collapsed[
  Cancer_v3$HPV_status_collapsed == "Non-tested"
] <- "Negative"

# Make it an ordered factor: Negative then Positive
Cancer_v3$HPV_status_collapsed <- factor(
  Cancer_v3$HPV_status_collapsed,
  levels = c("Negative", "Positive")
)

table(Cancer_v3$HPV_status_collapsed, useNA = "ifany")


hpv_cols <- c(
  "Negative" = "#4D9221",  # green
  "Positive" = "#C51B7D"   # magenta/red, matches tumor-ish color
)

############################################################
## 0. Choose reduction & extract UMAP limits (Cancer_v3)
############################################################

# Use the malignant-specific UMAP
umap_cancer <- "umap.cca_Cancer_v3"

umap_cancer_df <- Embeddings(Cancer_v3, umap_cancer)

xlim_cancer <- range(umap_cancer_df[, 1])
ylim_cancer <- range(umap_cancer_df[, 2])

# add 3% padding so points + stroke aren't clipped
pad_limits <- function(lim, frac = 0.03) {
  d <- diff(lim)
  c(lim[1] - d * frac, lim[2] + d * frac)
}

xlim_cancer_pad <- pad_limits(xlim_cancer, 0.03)
ylim_cancer_pad <- pad_limits(ylim_cancer, 0.03)


############################################################
## 1. Define cluster colors (your Okabe–Ito palette)
############################################################

# If you kept original labels: "Cancer_C0", ..., "Cancer_C4"
cluster_levels <- levels(factor(Cancer_v3$Unsupervised_Clusters))

cluster_colors <- setNames(
  c("#0072B2", "#009E73", "#CC79A7", "#F0E442", "#D55E00"),
  cluster_levels
)


# If you instead renamed to "Cluster 0"..."Cluster 4", then:
# cluster_levels <- paste("Cluster", 0:4)
# names(cluster_colors) <- cluster_levels

############################################################
## 2. Reuse your UMAP theme (no axes, legend ON)
############################################################

# You already defined this for RM_ST; reusing it is fine.
theme_umap_cancer <- function() {
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

############################################################
## 3. Cancer-specific UMAP function (mirrors RM_ST umap_png)
############################################################

cancer_umap_png <- function(group_by, colors, title_text) {
  
  DimPlot(
    Cancer_v3,
    reduction = umap_cancer,
    group.by  = group_by,
    pt.size   = 0.8,    # maybe tone this down a bit; 2.5 is huge
    raster    = F
  ) +
    scale_color_manual(values = colors, drop = FALSE) +
    coord_cartesian(
      xlim   = xlim_cancer_pad,
      ylim   = ylim_cancer_pad,
      expand = FALSE
    ) +
    theme_umap_cancer() +
    labs(
      x     = NULL,
      y     = NULL,
      color = title_text,
      title = title_text
    )
}



save_umap_with_legend <- function(group_by,
                                  colors,
                                  title_text,
                                  file_prefix,
                                  width  = 4,
                                  height = 3,
                                  dpi    = 900,
                                  legend_width = 2) {
  # 1. Full plot with legend
  p_full <- cancer_umap_png(
    group_by   = group_by,
    colors     = colors,
    title_text = title_text
  )
  
  # 2. UMAP panel only (no legend)
  p_nolegend <- p_full + theme(legend.position = "none")
  
  ggsave(
    paste0(file_prefix, ".png"),
    p_nolegend,
    width  = width,
    height = height,
    units  = "in",
    dpi    = dpi
  )
  
  # 3. Legend only
  leg <- cowplot::get_legend(
    p_full + theme(
      legend.position      = "right",
      legend.justification = "top"
    )
  )
  
  legend_plot <- cowplot::ggdraw(leg)
  
  ggsave(
    paste0(file_prefix, "_legend.png"),
    legend_plot,
    width  = legend_width,
    height = height,
    units  = "in",
    dpi    = dpi
  )
}
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Fig2abc/")  # removed (publication): use out_dir + file.path()
save_umap_with_legend(
  group_by   = "Unsupervised_Clusters",
  colors     = cluster_colors,  # your Okabe–Ito palette
  title_text = "Malignant Clusters (Unsupervised)",
  file_prefix = "UMAP_Cancer_UnsupervisedClusters"
)

FeaturePlot(Cancer_v3, features = "VIM", reduction = "umap.cca_Cancer_v3", pt.size = 1, min.cutoff = 1 )
FeaturePlot(Cancer_v3, features = "HIF1A", reduction = "umap.cca_Cancer_v3", pt.size = 1, max.cutoff = 2, min.cutoff = 0.5 )

FeaturePlot(Cancer_v3, features = "AUC_Hypoxia", reduction = "umap.cca_Cancer_v3", pt.size = 1, max.cutoff = 2, min.cutoff = 0.5 )

FeaturePlot(Cancer_v3, features = "AUC_Hypoxia", reduction = "umap.cca_Cancer_v3", pt.size = 1, max.cutoff = 2, min.cutoff = 0.5 )

FeaturePlot(Cancer_v3, features = "AUC_Hypoxia", reduction = "umap.cca_Cancer_v3", pt.size = 1, max.cutoff = 2, min.cutoff = 0.5 )

FeaturePlot(Cancer_v3, features = "AUC_EMT", reduction = "umap.cca_Cancer_v3", pt.size = 1, max.cutoff = 2, min.cutoff = 0.5 )

save_umap_with_legend(
  group_by   = "Sample_ID_renamed",
  colors     = sample_cols,   # the shuffled Dynamic palette
  title_text = "Sample ID",
  file_prefix = "UMAP_Cancer_SampleID"
)


save_umap_with_legend(
  group_by   = "HPV_status_collapsed",
  colors     = hpv_cols,  # e.g. Negative=green, Positive=magenta
  title_text = "HPV Status",
  file_prefix = "UMAP_Cancer_HPVstatus"
)



table(Cancer_v3$Cancer_Unsupervised_Clustering_v3_renamed)

DefaultAssay(Cancer_v3) <- "SCT"

Idents(Cancer_v3) <- "Unsupervised_Clusters"
Cancer_v3 <- BuildClusterTree(object = Cancer_v3, dims = 1:15,
                              reorder = FALSE,
                              reorder.numeric = FALSE)


tree <- Cancer_v3@tools$BuildClusterTree  # hclust
phy  <- as.phylo(tree)

phy$tip.label


# cluster levels from the tree (now aligned with your metadata)
cluster_levels <- phy$tip.label   # "Cluster 0"..."Cluster 4"

tree_dat <- as_tibble(phy)   # treeio::as_tibble.phylo

get_leaf_descendants <- function(node, phy) {
  phangorn::Descendants(phy, node, type = "tips")[[1]]
}

tree_dat$cluster <- sapply(tree_dat$node, function(node) {
  tip_ids <- get_leaf_descendants(node, phy)
  if (length(tip_ids) == 0) return(NA_character_)
  
  tip_labels <- phy$tip.label[tip_ids]
  tip_labels <- tip_labels[tip_labels %in% cluster_levels]
  if (length(tip_labels) == 0) return(NA_character_)
  
  tbl <- table(tip_labels)
  if (length(tbl) == 1) {
    return(names(tbl)[1])
  } else {
    return(names(tbl)[which.max(tbl)])  # majority cluster
  }
})

tree_dat$cluster <- factor(tree_dat$cluster, levels = cluster_levels)
tree_dat$cluster


# pretty labels = same as cluster_levels here
pretty_label_map <- setNames(cluster_levels, cluster_levels)
tree_dat$pretty_label <- pretty_label_map[tree_dat$label]

p <- ggtree(phy, aes(color = cluster)) %<+% tree_dat +
  geom_tree(size = 1.5) +
  ggtree::geom_point2(aes(subset = !is.na(cluster)), size = 4) +
  geom_tiplab(
    aes(label = pretty_label, color = cluster),
    offset   = 0.5,
    size     = 4,
    fontface = "bold"
  ) +
  scale_color_manual(
    values   = cluster_colors,
    na.value = "grey80",
    name     = "Malignant cluster"
  ) +
  scale_y_reverse() +          # ← ★★ FLIP TREE ON X-AXIS ★★
  coord_cartesian(clip = "off") +
  theme_tree() +
  theme(
    plot.margin     = unit(c(0.2, 2.5, 0.2, 0.2), "cm"),
    legend.position = "none"     # ← REMOVE LEGEND
  )

print(p)


ggsave(filename = file.path(out_dir, "Tree_Cancer_Unsupervised_clean.png"),
  p,           # your tree object from earlier
  width  = 4,
  height = 3,
  units  = "in",
  dpi    = 900
)



DefaultAssay(Cancer_v3) <- "SCT"   # or "RNA", wherever CLDN4/SPRR1B/etc live

p_neb <- SCpubr::do_NebulosaPlot(
  sample    = Cancer_v3,
  features  = c("CLDN4","SPRR1B","LAMC2","ITGA5"),
  reduction = "umap.cca_Cancer_v3",
  use_viridis          = FALSE,
  sequential.palette   = "YlGnBu",
  sequential.direction = 1,
  legend.position      = "bottom",
  legend.length        = 6,
  legend.width         = 0.4,
  plot.axes            = FALSE,
  font.size            = 12
)

ggsave(filename = file.path(out_dir, "Nebulosa_Markers_CLDN4_SPRR1B_LAMC2_ITGA5.png"),
  p_neb,
  width  = 7.5,
  height = 7,
  units  = "in",
  dpi    = 900
)



#DEGs -> Identifying Clusters: 

#DEGs -> Identifying Clusters: 
DEGs_Spatial_Cancer_filtered <- read.csv("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part2/2.3 - Visualization/DEGs_Cancer_V3.csv")

unique(DEGs_Spatial_Cancer_filtered$cluster)

# Replace 'cluster' with your actual column name if needed
clusters <- unique(DEGs_Spatial_Cancer_filtered$cluster)

DEGs_Spatial_Cancer_filtered <- DEGs_Spatial_Cancer_filtered %>%
  mutate(
    cluster = gsub("Cancer_C", "cluster ", cluster),
    cluster = ifelse(cluster == "cluster 0", "cluster 5", cluster)
  )

unique(DEGs_Spatial_Cancer_filtered$cluster)

# "Cluster 3" "Cluster 2" "Cluster 1" "Cluster 4" "Cluster 0"

DEGs_Spatial_Cancer_filtered$cluster <- factor(
  DEGs_Spatial_Cancer_filtered$cluster,
  levels = paste("cluster", 1:5)
)



DEGs_Spatial_Cancer_filtered <- DEGs_Spatial_Cancer_filtered %>%
  mutate(cluster = gsub("Cancer_C", "cluster ", cluster)) %>%
  mutate(cluster = factor(cluster, levels = paste("cluster", 1:5)))


# Get top 15 most upregulated genes per cluster

# Get top 15 upregulated genes per cluster
top15_up_per_cluster <- DEGs_Spatial_Cancer_filtered %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

# Extract gene names per cluster as a named list
top15_genes_list <- top15_up_per_cluster %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene)) %>%
  deframe()


# 1. Define your gene groups (ordered lists)
gene_groups <- top15_genes_list

# 2. Your color mapping per cluster
cancer_cluster_colors <- cluster_colors

# 3. Combine all genes to a unique gene list
gene_list <- unique(unlist(gene_groups))

# 4. Set assay and scale data for selected genes
DefaultAssay(Cancer_v3) <- "Spatial"
Cancer_v3 <- ScaleData(Cancer_v3, features = gene_list, verbose = FALSE)

# 5. Calculate average expression per group
df <- AverageExpression(
  Cancer_v3, assays = "Spatial", features = gene_list,
  group.by = "Unsupervised_Clusters", slot = "scale.data"
)$Spatial

df_t <- t(df)

# 6. Prepare and filter gene groups to keep only genes present in data
group_labels <- names(gene_groups)
gene_groups <- lapply(gene_groups, function(genes) genes[genes %in% colnames(df_t)])
non_empty_groups <- group_labels[sapply(gene_groups, length) > 0]

# 7. Map each gene to its group
gene2group <- c()
for (grp in names(gene_groups)) {
  g <- gene_groups[[grp]]
  gene2group <- c(gene2group, setNames(rep(grp, length(g)), g))
}

# 8. Define your consistent, desired gene order (grouped by cluster, original order)
desired_gene_order <- unlist(gene_groups[non_empty_groups])
desired_gene_order <- desired_gene_order[desired_gene_order %in% colnames(df_t)]

# 9. Reorder columns of your matrix and annotation accordingly
df_t <- df_t[, desired_gene_order, drop = FALSE]
col_groups <- gene2group[desired_gene_order]
annot_col <- data.frame(
  Top10Genes = factor(col_groups, levels = non_empty_groups)
)
rownames(annot_col) <- desired_gene_order

# 10. Calculate column gaps (for visual group separation)
gaps_col <- cumsum(sapply(non_empty_groups, function(g) sum(col_groups == g)))
if(length(gaps_col) > 1) gaps_col <- gaps_col[-length(gaps_col)]
if(length(gaps_col) == 0 || all(is.na(gaps_col))) gaps_col <- NULL

# 11. Define annotation colors
annotation_colors <- list(
  Top10Genes = cancer_cluster_colors[non_empty_groups]
)

colnames(annot_col) <- " "   # single space so it's basically invisible

annotation_colors <- list(
  " " = cancer_cluster_colors[non_empty_groups]
)

# 12. Plot heatmap

pheat <- pheatmap(
  df_t,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  annotation_col    = annot_col,
  annotation_colors = annotation_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  gaps_col          = gaps_col,
  angle_col         = 45,
  fontsize          = 15,
  fontsize_row      = 10,
  fontsize_col      = 10,
  annotation_legend = TRUE,
  legend            = TRUE,
  scale             = "none",
  color = colorRampPalette(c("#053061", "white", "#8b0000"))(500),
  breaks = seq(-1, 1, length.out = 501),
  
  
  # >>> THE KEY LINES <<<
  cellwidth  = 15,    # pixels (or "pt" units) — adjust to your needs
  cellheight = 15     # SAME VALUE → SQUARE BOXES
)


# (Optional) View heatmap
pheat


png(file.path(out_dir, "Cancer_Top15_Heatmap.png"),
    width  = 12000,   # 12 in × 200 dpi
    height = 1500,    # 4 in × 200 dpi
    res    = 600)

grid::grid.newpage()
grid::grid.draw(pheat$gtable)

dev.off()








###################################################
## 2. Top 15 upregulated genes per cluster
###################################################

top15_up_per_cluster <- DEGs_Spatial_Cancer_filtered %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15)

# Named list: "Cluster 0" -> c(gene1, gene2, ...)
gene_groups <- top15_up_per_cluster %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  deframe()

# sanity check
str(gene_groups)

###################################################
## 3. Color mapping per cluster
###################################################

cancer_cluster_colors <- cluster_colors     # from your earlier setup
# names(cancer_cluster_colors) should be "Cluster 0"..."Cluster 4"

###################################################
## 4. Build gene list & AverageExpression
###################################################

gene_list <- unique(unlist(gene_groups))

DefaultAssay(Cancer_v3) <- "Spatial"
Cancer_v3 <- ScaleData(Cancer_v3, features = gene_list, verbose = FALSE)

df <- AverageExpression(
  Cancer_v3,
  assays   = "Spatial",
  features = gene_list,
  group.by = "Unsupervised_Clusters",
  slot     = "scale.data"
)$Spatial
# df: rows = genes, cols = clusters
dim(df)
head(rownames(df))
colnames(df)    # should be "Cluster 0"... "Cluster 4"

###################################################
## 5. Filter gene_groups to genes present in df
###################################################

group_labels <- names(gene_groups)

gene_groups <- lapply(
  gene_groups,
  function(genes) genes[genes %in% rownames(df)]
)

non_empty_groups <- group_labels[sapply(gene_groups, length) > 0]
gene_groups   # should NOT be empty now

###################################################
## 6. Gene order + mapping gene -> cluster group
###################################################

gene2group <- c()
for (grp in names(gene_groups)) {
  g <- gene_groups[[grp]]
  gene2group <- c(gene2group, setNames(rep(grp, length(g)), g))
}

desired_gene_order <- unlist(gene_groups[non_empty_groups])
desired_gene_order <- desired_gene_order[desired_gene_order %in% rownames(df)]

# Reorder df rows: genes (rows) in desired order
df_flipped <- df[desired_gene_order, , drop = FALSE]

###################################################
## 7. Row annotation (left bar = gene’s cluster)
###################################################

row_groups <- gene2group[rownames(df_flipped)]

annot_row <- data.frame(
  GeneGroup = factor(row_groups, levels = non_empty_groups)
)
rownames(annot_row) <- rownames(df_flipped)

# If you want to hide the annotation title text:
colnames(annot_row) <- " "

annotation_colors <- list(
  " " = cancer_cluster_colors[non_empty_groups]
)

###################################################
## 8. Row gaps between cluster groups
###################################################

gaps_row <- cumsum(sapply(non_empty_groups, function(g) sum(row_groups == g)))
if (length(gaps_row) > 1) gaps_row <- gaps_row[-length(gaps_row)]
if (length(gaps_row) == 0 || all(is.na(gaps_row))) gaps_row <- NULL

###################################################
## 9. Flipped heatmap (genes as rows, clusters as cols)
###################################################

pheat <- pheatmap(
  df_flipped,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  annotation_row    = annot_row,
  annotation_colors = annotation_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  gaps_row          = gaps_row,
  angle_col         = 45,
  fontsize          = 15,
  fontsize_row      = 15,
  fontsize_col      = 15,
  annotation_legend = F,
  legend            = T,
  scale             = "none",
  color             = colorRampPalette(c("#053061", "white", "#8b0000"))(500),
  breaks            = seq(-1, 1, length.out = 501),
  cellwidth         = 18,
  cellheight        = 18,
  border_color      = "grey40",
)

png(file.path(out_dir, "Cancer_Top15_Heatmap_fliped.png"),
    width  = 3000,   # 12 in × 200 dpi
    height = 12000,    # 4 in × 200 dpi
    res    = 600)

grid::grid.newpage()
grid::grid.draw(pheat$gtable)

dev.off()


n_tiles <- 15
cellheight_pt <- 18
height_in <- n_tiles * cellheight_pt / 72  # 15 * 18 / 72 = 3.75 in





##AUC

# Read the Excel file
##These DEGs are obtained using DGE n1 vs n3 in Bose Lab



DEGs <- read_excel("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part2/2.0 - Malignant Spot Characterization/DEGs_BoseLabSignature.xlsx")

core <- subset(DEGs, group %in% "n1")
edge <- subset(DEGs, group %in% "n3")

core_genes <- core$gene #362
edge_genes <- edge$gene #120

intersected_genes_TC <- intersect(core_genes, rownames(Cancer_v3)) #349
intersected_genes_LE <- intersect(edge_genes, rownames(Cancer_v3)) #115


# HNSCC Cancer Cell State Markers
Differentiated_markers <- c("IVL", "SPRR1B", "SPRR2A", "SPRR2D", "TGM1", "TGM3", "KRT1", "KRT10", "LOR", "FLG")
Proliferative_markers <- c("MKI67", "PCNA", "TOP2A", "MCM2", "MCM4", "AURKA", "BIRC5")
EMT_markers <- c("VIM", "FN1", "SNAI2", "TWIST1", "ZEB1", "CDH2", "COL1A1", "COL1A2", "COL3A1", "COL5A1")
Hypoxic_markers <- c("CA9", "VEGFA","LDHA", "HIF1A", "ENO2", "NDRG1")


gene_sets <- list(
  Tumor_Core     = intersected_genes_TC,
  Leading_Edge   = intersected_genes_LE,
  Differentiation = Differentiated_markers,
  Proliferation   = Proliferative_markers,
  EMT            = EMT_markers,
  Hypoxia        = Hypoxic_markers
)

# Step 1: Get sparse raw counts from RNA assay
counts <- GetAssayData(Cancer_v3, assay = "Spatial", slot = "counts")

# Step 2: Compute column-wise (cell-wise) sum efficiently
cell_sums <- Matrix::colSums(counts)

# Step 3: Normalize to CPM manually (sparse-aware)
counts_cpm <- t(t(counts) / cell_sums * 10000)  # Still sparse!

# Step 4: Transpose to cells-as-rows for AUCell
exprMatrix <- as((counts_cpm), "dgCMatrix")    # now: cells x genes

# exprMatrix: cells x genes, already created earlier
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
cells_AUC      <- AUCell_calcAUC(gene_sets, cells_rankings)


# Extract AUC matrix: gene sets × cells
aucs   <- getAUC(cells_AUC)
aucs_t <- data.frame(t(aucs))  # cells × gene sets

# Attach cluster info
aucs_t$Unsupervised_Clusters <- Cancer_v3$Unsupervised_Clusters[match(rownames(aucs_t),
                                                                      colnames(Cancer_v3))]
# After:
# aucs_t <- data.frame(t(aucs))
# aucs_t$Unsupervised_Clusters <- Cancer_v3$Unsupervised_Clusters[match(rownames(aucs_t), colnames(Cancer_v3))]

# Names of your gene sets
sig_names <- names(gene_sets)
sig_names
# should be: "Tumor_Core" "Leading_Edge" "Differentiation" "Proliferation" "EMT" "Hypoxia"

# Find which of these are present as columns
old_names <- intersect(colnames(aucs_t), sig_names)
old_names
# should give the same 6 names above

# Rename those columns to add AUC_ prefix
colnames(aucs_t)[match(old_names, colnames(aucs_t))] <- paste0("AUC_", old_names)

# sanity check
colnames(aucs_t)[1:10]


aucs_meta <- aucs_t[, grepl("^AUC_", colnames(aucs_t))]
dim(aucs_meta)
# should be 9819 x 6 (rows = cells, cols = signatures)

Cancer_v3 <- AddMetaData(Cancer_v3, metadata = aucs_meta)

colnames(Cancer_v3@meta.data)[grepl("^AUC_", colnames(Cancer_v3@meta.data))]



# Order clusters as you like; adjust if you have more than 0–4
cluster_order <- paste("cluster", 1:5)


aucs_summary <- Cancer_v3@meta.data %>%
  dplyr::select(
    Unsupervised_Clusters,
    AUC_Tumor_Core,
    AUC_Leading_Edge,
    AUC_Differentiation,
    AUC_Proliferation,
    AUC_EMT,
    AUC_Hypoxia
  ) %>%
  group_by(Unsupervised_Clusters) %>%
  summarise(
    across(starts_with("AUC_"), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Unsupervised_Clusters = factor(Unsupervised_Clusters,
                                   levels = cluster_order)
  ) %>%
  filter(!is.na(Unsupervised_Clusters)) %>%
  arrange(Unsupervised_Clusters)

aucs_summary


aucs_summary_plot <- aucs_summary %>%
  rename_with(~ sub("^AUC_", "", .x), starts_with("AUC_"))
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Fig2abc//")  # removed (publication): use out_dir + file.path()

signature_colors <- c(
  Tumor_Core      = "#4DBBD5FF",
  Leading_Edge    = "#E64B35FF",
  Differentiation = "darkred",
  Proliferation   = "#117A65",
  EMT             = "#9966CC",
  Hypoxia         = "#1B4F72"
)

signatures <- names(signature_colors)
cluster_order <- paste("cluster", 1:5)

aucs_summary_plot$Unsupervised_Clusters <- factor(
  aucs_summary_plot$Unsupervised_Clusters,
  levels = rev(cluster_order)   # <-- reverse so Cluster 0 ends up at the top
)

# Pretty (publication) display names for plot titles
signature_pretty <- c(
  Tumor_Core      = "Tumor core gene signature",
  Leading_Edge    = "Leading edge gene signature",
  Differentiation = "Differentiation gene signature",
  Proliferation   = "Proliferation gene signature",
  EMT             = "EMT gene signature",
  Hypoxia         = "Hypoxia gene signature"
)

make_auc_barplot <- function(df, signature, color) {
  
  plot_title <- if (!is.null(signature_pretty[[signature]])) signature_pretty[[signature]] else signature
  
  ggplot(df, aes(x = Unsupervised_Clusters, y = .data[[signature]])) +
    geom_bar(stat = "identity", width = 0.7, fill = color) +
    geom_hline(
      yintercept = median(df[[signature]], na.rm = TRUE),
      linetype   = "dotted",
      color      = "black",
      linewidth  = 0.5
    ) +
    coord_flip() +
    ylab("Mean AUC") +
    xlab(NULL) +
    ggtitle(plot_title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(size = 12, face = "bold", hjust = 0),
      axis.text.y  = element_text(size = 10),
      axis.text.x  = element_text(size = 9),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
    )
}


auc_barplots <- lapply(signatures, function(sig) {
  make_auc_barplot(aucs_summary_plot, sig, signature_colors[sig])
})
names(auc_barplots) <- signatures

# Arrange 3x2 for a panel
p_bar <- wrap_plots(auc_barplots, ncol = 3)


ggsave(filename = file.path(out_dir, "AUC_cluster_barplots_6signatures.png"),
       plot = p_bar, width = 10, height = 6, units = "in", dpi = 300)



df_meta_long <- Cancer_v3@meta.data %>%
  dplyr::select(
    Unsupervised_Clusters,
    AUC_Tumor_Core,
    AUC_Leading_Edge,
    AUC_Differentiation,
    AUC_Proliferation,
    AUC_EMT,
    AUC_Hypoxia
  ) %>%
  mutate(
    Unsupervised_Clusters = factor(Unsupervised_Clusters, levels = cluster_order)
  ) %>%
  pivot_longer(
    cols = starts_with("AUC_"),
    names_to  = "Signature",
    values_to = "AUC"
  ) %>%
  mutate(
    Signature = sub("^AUC_", "", Signature),
    Signature = factor(Signature,
                       levels = c("Tumor_Core", "Leading_Edge",
                                  "Differentiation", "Proliferation",
                                  "EMT", "Hypoxia"))
  ) %>%
  filter(!is.na(Unsupervised_Clusters), !is.na(AUC))

p_violin_box <- ggplot(df_meta_long,
                       aes(x = Unsupervised_Clusters,
                           y = AUC,
                           fill = Unsupervised_Clusters)) +
  
  geom_violin(
    trim  = TRUE,
    scale = "width",
    width = 0.9,
    color = NA,
    alpha = 0.85
  ) +
  geom_boxplot(
    width  = 0.25,
    outlier.size = 0.3,
    outlier.alpha = 0.4,
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = cluster_colors) +
  
  facet_wrap(~ Signature, ncol = 3, scales = "free_y") +
  xlab(NULL) +
  ylab("AUCell score") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(size = 11, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position  = "none"
  )

p_violin_box
ggsave(filename = file.path(out_dir, "AUC_violin_box_by_cluster_6signatures.png"),
       plot = p_violin_box, width = 10, height = 6, units = "in", dpi = 300)






cluster_order <- paste("Cluster", 0:4)

df_meta_long <- Cancer_v3@meta.data %>%
  dplyr::select(
    Unsupervised_Clusters,
    AUC_Tumor_Core,
    AUC_Leading_Edge,
    AUC_Differentiation,
    AUC_Proliferation,
    AUC_EMT,
    AUC_Hypoxia
  ) %>%
  mutate(
    Unsupervised_Clusters = factor(Unsupervised_Clusters, levels = cluster_order)
  ) %>%
  pivot_longer(
    cols      = starts_with("AUC_"),
    names_to  = "Signature",
    values_to = "AUC"
  ) %>%
  mutate(
    Signature = sub("^AUC_", "", Signature),
    Signature = factor(
      Signature,
      levels = c("Tumor_Core", "Leading_Edge",
                 "Differentiation", "Proliferation",
                 "EMT", "Hypoxia")
    )
  ) %>%
  filter(!is.na(Unsupervised_Clusters), !is.na(AUC))




cluster_order <- paste("Cluster", 0:4)

df_meta_long$Unsupervised_Clusters <- factor(
  df_meta_long$Unsupervised_Clusters,
  levels = rev(cluster_order)     # <-- THIS FIXES IT
)

cluster_colors <- c(
  "Cluster 0" = "#D55E00",
  "Cluster 1" = "#0072B2",
  "Cluster 2" = "#009E73",
  "Cluster 3" = "#CC79A7",
  "Cluster 4" = "#F0E442"
)

p_violin_horiz <- ggplot(
  df_meta_long,
  aes(x = AUC, y = Unsupervised_Clusters, fill = Unsupervised_Clusters)
) +
  geom_violin(trim = TRUE, scale = "width", width = 0.9,
              color = NA, alpha = 0.85) +
  geom_boxplot(width = 0.25, outlier.size = 0.3, outlier.alpha = 0.4,
               fill = "white", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cluster_colors) +
  facet_grid(Signature ~ ., scales = "free_x") +
  xlab("AUCell score") + ylab(NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(size = 11, face = "bold"),
    axis.text.x      = element_text(size = 9),
    axis.text.y      = element_text(size = 10),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position  = "none"
  )

p_violin_horiz



ggsave(filename = file.path(out_dir, "AUC_horizontal_violin_box_by_cluster_6signatures.png"),
  plot = p_violin_horiz,
  width = 5, height = 10, units = "in",
  dpi = 600
)