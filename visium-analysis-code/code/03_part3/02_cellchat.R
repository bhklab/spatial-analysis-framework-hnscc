#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(CellChat)
  library(patchwork)  # CellChat sometimes loads it internally; safe to keep
})

# ============================================================
# Part 3 / CellChat — Tier1, Tier2, Tier2b, Tier2c 

# Tier1  : malignant-only (Cancer object)
# Tier2  : global (RM_ST with relabeled/condensed groups)
# Tier2b : malignant + CAF
# Tier2c : malignant + immune
#
# Outputs:
#   - *.rds CellChat objects
#   - communication tables (CSV)
#   - net matrices (CSV)
#
# IMPORTANT: no visualization here (handled in Part 4).
# ============================================================

# -------------------------
# EDIT PATHS
# -------------------------
rm_st_rds   <- "part3/RM_ST.rds"
cancer_rds  <- "part3/Cancer_v3.RDS"   # malignant-only object used for Tier1
out_root    <- "part3/cellchat"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# GLOBAL PARAMETERS
# -------------------------
assay_use          <- "Spatial"
slot_use           <- "data"      # CellChat typically uses normalized "data"
species_db         <- "human"     # change to "mouse" if needed
min_cells_per_group <- 10         # filter out tiny groups for stability
n_workers          <- 1           # keep 1 for reproducibility unless you parallelize explicitly

# -------------------------
# GROUPING COLUMNS
# -------------------------
# Tier1 grouping column in Cancer object:
tier1_group_by_candidates <- c("Final_Malignant_SubCluster",
                               "Cancer_Unsupervised_Clustering_v3",
                               "Cancer_SubClusters")
# Tier2 uses a new relabel column in RM_ST:
tier2_group_by <- "Spatial_CT_relabel"

# -------------------------
# RELABEL MAP (Tier2)
#   You can TURN THIS OFF by setting use_relabel_map <- FALSE
#   If OFF: Tier2 will use RM_ST$Spatial_CT directly as Spatial_CT_relabel.
# -------------------------
use_relabel_map <- TRUE

# This map is based on your current script; edit freely.
# If you want to postpone biological naming to Part 4, set use_relabel_map <- FALSE.
rename_map <- c(
  # malignant (example; turn off if delaying naming)
  "Neuroendocrine/Stem-like Niche" = "cycling",
  "TC_Epithelial_Differentiated"   = "TC",
  "LE_Pro_Inflammatory"            = "neutrophil_inflamed",
  "LE_ECM_Remodeling"              = "fibrovascular",
  "Hypoxic/Angiogenic Niche"       = "pEMT",
  
  # CAF (example)
  "Hypoxic/Invasion-associated CAFs (transitional)" = "non-activated fibroblast",
  "mCAFs"                                           = "ecm_myCAF",
  "apCAFs"                                          = "IFNg-iCAF",
  "iCAFs"                                           = "IL-iCAF",
  "Myofibroblast-like CAFs"                         = "acto-myCAF",
  
  # TAM / myeloid (example)
  "IFN-TAM"                   = "IFN-TAM",
  "pro-inflammatory TAM"      = "pro-inflammatory TAM",
  "monocyte"                  = "monocyte",
  
  # immune / other (example)
  "T cells"                   = "T cells",
  "Dendritic cells"           = "Dendritic cells",
  "Endothelial cells"         = "Endothelial cells",
  "antibody-secreting B cell" = "antibody-secreting B cell"
)

# -------------------------
# Tier2b / Tier2c GROUP SETS (must match tier2_group_by labels)
#   If use_relabel_map <- FALSE, update these to match your original labels.
# -------------------------
malig_groups <- c("cycling", "TC", "neutrophil_inflamed", "fibrovascular", "pEMT")

caf_groups <- c(
  "non-activated fibroblast",
  "ecm_myCAF",
  "IFNg-iCAF",
  "IL-iCAF",
  "acto-myCAF"
)

immune_groups <- c(
  "IFN-TAM",
  "pro-inflammatory TAM",
  "monocyte",
  "T cells",
  "Dendritic cells",
  "Endothelial cells",
  "antibody-secreting B cell"
)

# ============================================================
# Helpers
# ============================================================

pick_group_by <- function(obj, candidates) {
  md <- obj@meta.data
  for (nm in candidates) {
    if (nm %in% colnames(md)) return(nm)
  }
  stop("None of the requested group.by columns exist: ", paste(candidates, collapse = ", "))
}

prep_seurat_for_cellchat <- function(obj, assay = "Spatial", slot = "data") {
  DefaultAssay(obj) <- assay
  
  # Drop images to reduce size and avoid serialization issues
  if (length(obj@images) > 0) obj@images <- list()
  
  # Ensure normalized data exists
  # In Seurat v5, "data" may be layer-dependent; JoinLayers helps.
  if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    suppressWarnings({
      obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
    })
  }
  
  mat <- tryCatch(GetAssayData(obj, assay = assay, slot = slot),
                  error = function(e) NULL)
  
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("No usable ", slot, " slot detected. Running NormalizeData() on assay=", assay)
    obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
  }
  
  obj
}

run_cellchat_core <- function(obj, group.by, out_dir,
                              assay = "Spatial", slot = "data",
                              min_cells = 10, species_db = "human") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare expression matrix
  DefaultAssay(obj) <- assay
  data.input <- GetAssayData(obj, assay = assay, slot = slot)
  
  # Metadata
  meta <- obj@meta.data
  if (!(group.by %in% colnames(meta))) stop("Missing group.by column: ", group.by)
  meta[[group.by]] <- as.character(meta[[group.by]])
  
  # Filter rare groups (stabilizes computeCommunProb)
  keep_groups <- names(which(table(meta[[group.by]]) >= min_cells))
  meta[[group.by]] <- ifelse(meta[[group.by]] %in% keep_groups, meta[[group.by]], NA_character_)
  keep_cells <- rownames(meta)[!is.na(meta[[group.by]])]
  if (length(keep_cells) < 50) stop("Too few cells/spots after filtering by min_cells. Lower min_cells or check metadata.")
  
  data.input <- data.input[, keep_cells, drop = FALSE]
  meta <- meta[keep_cells, , drop = FALSE]
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = group.by)
  
  # Database
  if (species_db == "human") {
    cellchat@DB <- CellChatDB.human
  } else if (species_db == "mouse") {
    cellchat@DB <- CellChatDB.mouse
  } else {
    stop("species_db must be 'human' or 'mouse'")
  }
  
  # Standard pipeline (no plots)
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # You can tweak these if needed; leaving defaults for reproducibility
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Save object
  saveRDS(cellchat, file.path(out_dir, "cellchat.rds"))
  
  # Export communications
  comm <- subsetCommunication(cellchat)
  write.csv(comm, file.path(out_dir, "communications_all.csv"), row.names = FALSE)
  
  # Export net matrices
  if (!is.null(cellchat@net$count)) {
    write.csv(cellchat@net$count, file.path(out_dir, "net_count_matrix.csv"))
  }
  if (!is.null(cellchat@net$weight)) {
    write.csv(cellchat@net$weight, file.path(out_dir, "net_weight_matrix.csv"))
  }
  
  # Also export pathway-level communications (netP)
  # Some CellChat versions store netP as list; we export what’s easily tabular
  commP <- tryCatch(subsetCommunication(cellchat, slot.name = "netP"),
                    error = function(e) NULL)
  if (!is.null(commP)) {
    write.csv(commP, file.path(out_dir, "communications_pathway_level.csv"), row.names = FALSE)
  }
  
  message("DONE: ", out_dir)
  invisible(cellchat)
}

# ============================================================
# Load objects
# ============================================================
RM_ST <- readRDS(rm_st_rds)
Cancer <- readRDS(cancer_rds)

RM_ST <- prep_seurat_for_cellchat(RM_ST, assay = assay_use, slot = slot_use)
Cancer <- prep_seurat_for_cellchat(Cancer, assay = assay_use, slot = slot_use)

# ============================================================
# Tier1: malignant-only CellChat
# ============================================================
tier1_group_by <- pick_group_by(Cancer, tier1_group_by_candidates)
message("Tier1 group.by = ", tier1_group_by)

run_cellchat_core(
  obj      = Cancer,
  group.by = tier1_group_by,
  out_dir  = file.path(out_root, "tier1_malignant_only"),
  assay    = assay_use,
  slot     = slot_use,
  min_cells = min_cells_per_group,
  species_db = species_db
)

# ============================================================
# Tier2: global relabel + CellChat
# ============================================================
if (!("Spatial_CT" %in% colnames(RM_ST@meta.data))) {
  stop("RM_ST missing Spatial_CT. Provide Spatial_CT or update to your label column.")
}

orig_ct <- as.character(RM_ST$Spatial_CT)
Spatial_CT_relabel <- orig_ct

if (use_relabel_map) {
  to_rename <- Spatial_CT_relabel %in% names(rename_map)
  Spatial_CT_relabel[to_rename] <- rename_map[Spatial_CT_relabel[to_rename]]
}

RM_ST[[tier2_group_by]] <- Spatial_CT_relabel

run_cellchat_core(
  obj      = RM_ST,
  group.by = tier2_group_by,
  out_dir  = file.path(out_root, "tier2_global"),
  assay    = assay_use,
  slot     = slot_use,
  min_cells = min_cells_per_group,
  species_db = species_db
)

# ============================================================
# Tier2b: malignant + CAF only
# ============================================================
Idents(RM_ST) <- tier2_group_by
keep_tier2b <- c(malig_groups, caf_groups)
RM_ST_caf <- subset(RM_ST, idents = keep_tier2b)

run_cellchat_core(
  obj      = RM_ST_caf,
  group.by = tier2_group_by,
  out_dir  = file.path(out_root, "tier2b_malignant_plus_caf"),
  assay    = assay_use,
  slot     = slot_use,
  min_cells = min_cells_per_group,
  species_db = species_db
)

# ============================================================
# Tier2c: malignant + immune only
# ============================================================
Idents(RM_ST) <- tier2_group_by
keep_tier2c <- c(malig_groups, immune_groups)
RM_ST_immune <- subset(RM_ST, idents = keep_tier2c)

run_cellchat_core(
  obj      = RM_ST_immune,
  group.by = tier2_group_by,
  out_dir  = file.path(out_root, "tier2c_malignant_plus_immune"),
  assay    = assay_use,
  slot     = slot_use,
  min_cells = min_cells_per_group,
  species_db = species_db
)

