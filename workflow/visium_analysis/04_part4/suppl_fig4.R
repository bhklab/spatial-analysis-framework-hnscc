# ============================================================
# Supplementary Figure 4
# Publication polished script (paths parameterized, no setwd)
# ============================================================


# ---- shared helpers ----
# Assumes this script lives in the same folder as 00_common.R.
# If not, change the path below (e.g., source("code/04_part4/00_common.R")).
source("00_common.R")
theme_set(theme_pub())

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(writexl)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------


# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "suppl_fig4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


###Re-making the KEGG/Hallmark
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part2/2.3 - Visualization/")  # removed (publication): use out_dir + file.path()
if (!exists("DEGs_Spatial_caf_filtered")) {
  DEGs_Spatial_caf_filtered <- read_csv("caf_deg_outputs/CAFs_DEGs_filtered.csv",
                                        show_col_types = FALSE)
}
if (!exists("DEGs_Spatial_tam_filtered")) {
  DEGs_Spatial_tam_filtered <- read_csv("tam_deg_outputs/TAM_DEGs__tTAMs__filtered.csv",
                                        show_col_types = FALSE)
}

stopifnot(all(c("cluster","gene") %in% colnames(DEGs_Spatial_caf_filtered)))
stopifnot(all(c("cluster","gene") %in% colnames(DEGs_Spatial_tam_filtered)))


caf_cluster_remap <- c(
  "Hypoxic/Invasion-associated CAFs (transitional)"    = "non-activated fibroblast",
  "mCAFs"                                              = "ecm-myCAF",
  "apCAFs"                                             = "IFNg-iCAF",
  "Hypoxic/Invasion-associated CAFs (high activation)" = "IL-iCAF",
  "myCAFs"                                             = "acto-myCAF",
  "Erythroid/Platelet-interacting CAFs"                = "Erythroid+Platelet interacting CAFs"
)



###############################################
## CAF-only ORA (KEGG + Hallmark) + Directional ORA
## With renamed CAF cluster labels
###############################################

## -------------------------
## 0) Libraries
## -------------------------

###############################################
## 1) Paths
###############################################

# Where the DEGs + original outputs live
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part2/2.3 - Visualization/")  # removed (publication): use out_dir + file.path()
# Where you want all CAF ORA outputs to go (Fig4ab)
caf_out_base     <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/CAF_ORA"
caf_dir_out      <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/CAF_ORA_directional"
dir.create(caf_out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(caf_dir_out,  recursive = TRUE, showWarnings = FALSE)

###############################################
## 2) Load CAF DEGs (non-directional ORA)
###############################################

if (!exists("DEGs_Spatial_caf_filtered")) {
  DEGs_Spatial_caf_filtered <- read_csv(
    "caf_deg_outputs/CAFs_DEGs_filtered.csv",
    show_col_types = FALSE
  )
}

stopifnot(all(c("cluster", "gene") %in% colnames(DEGs_Spatial_caf_filtered)))

###############################################
## 3) CAF cluster renaming (for ORA)
###############################################

caf_cluster_remap <- c(
  "Hypoxic/Invasion-associated CAFs (transitional)"    = "non-activated fibroblast",
  "mCAFs"                                              = "ecm-myCAF",
  "apCAFs"                                             = "IFNg-iCAF",
  "Hypoxic/Invasion-associated CAFs (high activation)" = "IL-iCAF",
  "myCAFs"                                             = "acto-myCAF",
  "Erythroid/Platelet-interacting CAFs"                = "Erythroid+Platelet interacting CAFs"
)

caf_levels <- c(
  "non-activated fibroblast",
  "ecm-myCAF",
  "IFNg-iCAF",
  "IL-iCAF",
  "acto-myCAF",
  "Erythroid+Platelet interacting CAFs"
)

DEGs_Spatial_caf_filtered$cluster_original <- DEGs_Spatial_caf_filtered$cluster

DEGs_Spatial_caf_filtered$cluster <- caf_cluster_remap[DEGs_Spatial_caf_filtered$cluster]
DEGs_Spatial_caf_filtered$cluster <- factor(DEGs_Spatial_caf_filtered$cluster,
                                            levels = caf_levels)

table(DEGs_Spatial_caf_filtered$cluster, useNA = "ifany")

###############################################
## 4) Hallmark gene sets + mapping helper
###############################################

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

map_symbols_to_entrez <- function(symbols) {
  if (length(symbols) == 0) return(character(0))
  mapped <- bitr(symbols,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db,
                 drop     = TRUE)
  unique(mapped$ENTREZID)
}

###############################################
## 5) Core ORA runner (KEGG + Hallmark) per cluster
###############################################

run_enrichment_for_table <- function(deg_tbl,
                                     out_dir,
                                     prefix = "CAF",
                                     kegg_min_size = 10,
                                     kegg_max_size = 500,
                                     hallmark_min_size = 10) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # split DEGs by cluster
  deg_by_cluster <- split(deg_tbl, deg_tbl$cluster)
  
  kegg_sheets <- list()
  hall_sheets <- list()
  
  for (cl in names(deg_by_cluster)) {
    dat <- deg_by_cluster[[cl]]
    
    symbols <- unique(dat$gene)
    entrez  <- map_symbols_to_entrez(symbols)
    
    if (length(entrez) < hallmark_min_size) next
    
    ## KEGG ORA
    ekegg <- tryCatch({
      enrichKEGG(gene          = entrez,
                 organism      = "hsa",
                 keyType       = "ncbi-geneid",
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.25,
                 minGSSize     = kegg_min_size,
                 maxGSSize     = kegg_max_size)
    }, error = function(e) NULL)
    
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      kegg_res <- as.data.frame(ekegg) %>%
        mutate(cluster = cl, .before = 1)
      kegg_sheets[[paste0("KEGG_", gsub("[^A-Za-z0-9]+","_", cl))]] <- kegg_res
      
      readr::write_csv(
        kegg_res,
        file.path(out_dir, paste0(prefix, "__", gsub("[/ ]","+ ", cl), "__KEGG.csv"))
      )
      
      # quick dotplot
      try({
        dp <- dotplot(ekegg, showCategory = 15) +
          ggtitle(paste(prefix, "KEGG:", cl))
        ggsave(
          file.path(out_dir, paste0(prefix, "__", gsub("[/ ]","+ ", cl), "__KEGG_dotplot.pdf")),
          dp, width = 7, height = 5
        )
      }, silent = TRUE)
    }
    
    ## Hallmark ORA (SYMBOL + TERM2GENE)
    eh <- tryCatch({
      enricher(gene          = symbols,
               TERM2GENE     = hallmark,
               pAdjustMethod = "BH")
    }, error = function(e) NULL)
    
    if (!is.null(eh) && nrow(as.data.frame(eh)) > 0) {
      hall_res <- as.data.frame(eh) %>%
        mutate(cluster = cl, .before = 1)
      hall_sheets[[paste0("HALL_", gsub("[^A-Za-z0-9]+","_", cl))]] <- hall_res
      
      readr::write_csv(
        hall_res,
        file.path(out_dir, paste0(prefix, "__", gsub("[/ ]","+ ", cl), "__HALLMARK.csv"))
      )
      
      try({
        dp2 <- dotplot(eh, showCategory = 15) +
          ggtitle(paste(prefix, "Hallmark:", cl))
        ggsave(
          file.path(out_dir, paste0(prefix, "__", gsub("[/ ]","+ ", cl), "__HALLMARK_dotplot.pdf")),
          dp2, width = 7, height = 5
        )
      }, silent = TRUE)
    }
  }
  
  ## Excel workbook
  wb_list <- list()
  if (length(kegg_sheets)) {
    combined_kegg <- dplyr::bind_rows(kegg_sheets, .id = "sheet_src")
    wb_list$KEGG__combined <- combined_kegg
    wb_list <- c(wb_list, kegg_sheets)
  }
  if (length(hall_sheets)) {
    combined_hall <- dplyr::bind_rows(hall_sheets, .id = "sheet_src")
    wb_list$HALLMARK__combined <- combined_hall
    wb_list <- c(wb_list, hall_sheets)
  }
  
  if (length(wb_list)) {
    writexl::write_xlsx(
      wb_list,
      file.path(out_dir, paste0(prefix, "__ORA_results.xlsx"))
    )
  }
  
  message("✓ ORA finished for ", prefix, " → ", out_dir)
}

###############################################
## 6) Run standard ORA for CAFs (renamed clusters)
###############################################

run_enrichment_for_table(
  deg_tbl = DEGs_Spatial_caf_filtered,
  out_dir = caf_out_base,
  prefix  = "CAF"
)

###############################################
## 7) Directional ORA (Up vs Down) for CAFs
###############################################

###############################################
## Helper functions for directional ORA + plotting
## (run once per session, ABOVE step 7)
###############################################


## ---------- Hallmark TERM2GENE ----------

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_t2g <- hallmark_df %>%
  dplyr::select(
    term = dplyr::any_of(c("gs_name", "gs_name")),
    gene = dplyr::any_of(c("gene_symbol", "human_gene_symbol"))
  ) %>%
  distinct()

## ---------- small helpers ----------

map_symbols_to_entrez <- function(symbols) {
  if (length(symbols) == 0) return(character(0))
  out <- bitr(unique(symbols),
              fromType = "SYMBOL",
              toType   = "ENTREZID",
              OrgDb    = org.Hs.eg.db,
              drop     = TRUE)
  unique(out$ENTREZID)
}

safe_name <- function(x) gsub("[/ ]","+ ", x)

## ==========================================================
##  run_directional_ora_top10
## ==========================================================
run_directional_ora_top10 <- function(deg_tbl,
                                      out_dir,
                                      prefix = c("CAF","TAM"),
                                      kegg_min_size = 10,
                                      kegg_max_size = 500,
                                      hall_min_size = 10) {
  prefix <- match.arg(prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  stopifnot(all(c("cluster","gene","avg_log2FC") %in% colnames(deg_tbl)))
  
  # Universe (SYMBOL) = all genes present in this DEG table
  universe_symbols <- unique(deg_tbl$gene)
  universe_entrez  <- map_symbols_to_entrez(universe_symbols)
  
  # Split DEGs by cluster label
  by_cl <- split(deg_tbl, deg_tbl$cluster)
  
  # Collect all top10 results for Excel/combined CSV
  top10_collect <- list()
  
  sheet_nm <- function(base) {
    paste0(substr(gsub("[^A-Za-z0-9]+","_", base), 1, 28))
  }
  
  for (cl in names(by_cl)) {
    dat <- by_cl[[cl]]
    if (!nrow(dat)) next
    
    # Directional splits
    up_genes   <- dat %>% filter(avg_log2FC > 0) %>% pull(gene) %>% unique()
    down_genes <- dat %>% filter(avg_log2FC < 0) %>% pull(gene) %>% unique()
    
    # ---- KEGG (ENTREZ) ----
    up_entrez   <- map_symbols_to_entrez(up_genes)
    down_entrez <- map_symbols_to_entrez(down_genes)
    
    kegg_up <- kegg_down <- NULL
    
    if (length(up_entrez) >= kegg_min_size) {
      kegg_up <- tryCatch(
        enrichKEGG(
          gene          = up_entrez,
          universe      = universe_entrez,
          organism      = "hsa",
          keyType       = "ncbi-geneid",
          pAdjustMethod = "BH",
          minGSSize     = kegg_min_size,
          maxGSSize     = kegg_max_size
        ),
        error = function(e) NULL
      )
    }
    
    if (length(down_entrez) >= kegg_min_size) {
      kegg_down <- tryCatch(
        enrichKEGG(
          gene          = down_entrez,
          universe      = universe_entrez,
          organism      = "hsa",
          keyType       = "ncbi-geneid",
          pAdjustMethod = "BH",
          minGSSize     = kegg_min_size,
          maxGSSize     = kegg_max_size
        ),
        error = function(e) NULL
      )
    }
    
    # ---- Hallmark (SYMBOL) ----
    hall_up <- hall_down <- NULL
    
    if (length(up_genes) >= hall_min_size) {
      hall_up <- tryCatch(
        enricher(
          gene          = up_genes,
          universe      = universe_symbols,
          TERM2GENE     = hallmark_t2g,
          pAdjustMethod = "BH"
        ),
        error = function(e) NULL
      )
    }
    
    if (length(down_genes) >= hall_min_size) {
      hall_down <- tryCatch(
        enricher(
          gene          = down_genes,
          universe      = universe_symbols,
          TERM2GENE     = hallmark_t2g,
          pAdjustMethod = "BH"
        ),
        error = function(e) NULL
      )
    }
    
    # Convert to data.frame, add direction + cluster + score
    dfize <- function(x, direction) {
      if (is.null(x)) return(NULL)
      as.data.frame(x) %>%
        mutate(
          direction = direction,
          cluster   = cl,
          .before   = 1
        ) %>%
        mutate(
          score = -log10(p.adjust),
          .after = p.adjust
        )
    }
    
    kegg_up_df   <- dfize(kegg_up,   "Activated (Up)")
    kegg_down_df <- dfize(kegg_down, "Inhibited (Down)")
    hall_up_df   <- dfize(hall_up,   "Activated (Up)")
    hall_down_df <- dfize(hall_down, "Inhibited (Down)")
    
    # Top 10 per direction by adjusted p-value
    top10_kegg <- bind_rows(kegg_up_df, kegg_down_df) %>%
      arrange(direction, p.adjust) %>%
      group_by(direction) %>%
      slice_head(n = 10) %>%
      ungroup()
    
    top10_hall <- bind_rows(hall_up_df, hall_down_df) %>%
      arrange(direction, p.adjust) %>%
      group_by(direction) %>%
      slice_head(n = 10) %>%
      ungroup()
    
    # Save per-cluster CSVs
    if (nrow(top10_kegg)) {
      write_csv(
        top10_kegg,
        file.path(out_dir,
                  paste0(prefix, "__", safe_name(cl),
                         "__KEGG_top10_directional.csv"))
      )
    }
    
    if (nrow(top10_hall)) {
      write_csv(
        top10_hall,
        file.path(out_dir,
                  paste0(prefix, "__", safe_name(cl),
                         "__HALLMARK_top10_directional.csv"))
      )
    }
    
    # Collect for combined Excel
    if (nrow(top10_kegg)) {
      top10_collect[[paste0("KEGG__", sheet_nm(cl))]] <- top10_kegg
    }
    if (nrow(top10_hall)) {
      top10_collect[[paste0("HALL__", sheet_nm(cl))]] <- top10_hall
    }
  }
  
  # Combined CSV + Excel workbook
  if (length(top10_collect)) {
    combined <- bind_rows(top10_collect, .id = "sheet_src")
    write_csv(
      combined,
      file.path(out_dir,
                paste0(prefix, "__directional_ORA_top10__combined.csv"))
    )
    
    wb_list <- list()
    kegg_all <- combined %>%
      filter(grepl("^KEGG__", sheet_src))
    hall_all <- combined %>%
      filter(grepl("^HALL__", sheet_src))
    
    if (nrow(kegg_all)) wb_list$KEGG__all_clusters <- kegg_all
    if (nrow(hall_all)) wb_list$HALLMARK__all_clusters <- hall_all
    
    wb_list <- c(wb_list, top10_collect)
    
    write_xlsx(
      wb_list,
      file.path(out_dir,
                paste0(prefix, "__directional_ORA_top10.xlsx"))
    )
  }
  
  message("✓ Directional ORA top10 done for ", prefix, " → ", out_dir)
}

## ==========================================================
## Plotting helper: plot_dir_ora_top10
## ==========================================================

.read_top10 <- function(dir_path, prefix, db = c("KEGG","HALLMARK")) {
  db <- match.arg(db)
  combo <- file.path(dir_path,
                     paste0(prefix, "__directional_ORA_top10__combined.csv"))
  
  .filter_db <- function(df, db) {
    if (!"sheet_src" %in% names(df)) return(df)
    if (db == "KEGG") {
      dplyr::filter(df, grepl("^KEGG__", .data$sheet_src))
    } else {
      dplyr::filter(df,
                    grepl("^HALL(__|MARK__)", .data$sheet_src) |
                      grepl("^HALL__", .data$sheet_src))
    }
  }
  
  if (file.exists(combo)) {
    res <- readr::read_csv(combo, show_col_types = FALSE) |> .filter_db(db)
    if (nrow(res) > 0) return(res)
  }
  
  pat <- if (db == "KEGG") {
    paste0("^", prefix, "__.*KEGG_top10_directional\\.csv$")
  } else {
    paste0("^", prefix, "__.*(HALLMARK|HALL)_top10_directional\\.csv$")
  }
  files <- list.files(dir_path, pattern = pat, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))
}

.plot_cluster <- function(df_cluster, title_prefix, outfile) {
  if (is.null(df_cluster) || nrow(df_cluster) == 0) return(invisible(NULL))
  
  dfp <- df_cluster %>%
    dplyr::select(cluster, direction, Description, p.adjust, score) %>%
    dplyr::mutate(Description = as.character(Description)) %>%
    dplyr::group_by(direction) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::mutate(
      Description_fac = forcats::fct_reorder(Description, score)
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot(dfp, aes(x = Description_fac, y = score)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ direction, ncol = 1, scales = "free_y") +
    labs(
      title = paste0(title_prefix, ": ", unique(dfp$cluster)),
      x     = NULL,
      y     = expression(-log[10]("adj. p"))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold")
    )
  
  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

plot_dir_ora_top10 <- function(dir_path, prefix) {
  message("Plotting directional ORA top10 for: ", prefix)
  
  # KEGG
  kegg_df <- .read_top10(dir_path, prefix, "KEGG")
  if (!is.null(kegg_df) && nrow(kegg_df)) {
    out_kegg <- file.path(dir_path, "plots_kegg")
    dir.create(out_kegg, showWarnings = FALSE, recursive = TRUE)
    for (cl in sort(unique(kegg_df$cluster))) {
      this <- dplyr::filter(kegg_df, cluster == cl)
      .plot_cluster(
        this,
        paste0(prefix, " KEGG"),
        file.path(out_kegg,
                  paste0(prefix, "__", safe_name(cl),
                         "__KEGG_top10_directional.png"))
      )
    }
  } else {
    message("No KEGG top10 data found in: ", dir_path)
  }
  
  # Hallmark
  hall_df <- .read_top10(dir_path, prefix, "HALLMARK")
  if (!is.null(hall_df) && nrow(hall_df)) {
    out_hall <- file.path(dir_path, "plots_hallmark")
    dir.create(out_hall, showWarnings = FALSE, recursive = TRUE)
    for (cl in sort(unique(hall_df$cluster))) {
      this <- dplyr::filter(hall_df, cluster == cl)
      .plot_cluster(
        this,
        paste0(prefix, " Hallmark"),
        file.path(out_hall,
                  paste0(prefix, "__", safe_name(cl),
                         "__HALLMARK_top10_directional.png"))
      )
    }
  } else {
    message("No Hallmark top10 data found in: ", dir_path)
  }
}

###############################################
## 7) Directional ORA (Up vs Down) for CAFs
##    (KEGG + Hallmark, using renamed CAF clusters)
###############################################


## Assumes you already ran:
##   - setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part2/2.3 - Visualization/")
##   - DEGs_Spatial_caf_filtered <- read_csv("caf_deg_outputs/CAFs_DEGs_filtered.csv", show_col_types = FALSE)
##   - caf_cluster_remap <- c(
##       "Hypoxic/Invasion-associated CAFs (transitional)"    = "non-activated fibroblast",
##       "mCAFs"                                              = "ecm-myCAF",
##       "apCAFs"                                             = "IFNg-iCAF",
##       "Hypoxic/Invasion-associated CAFs (high activation)" = "IL-iCAF",
##       "myCAFs"                                             = "acto-myCAF",
##       "Erythroid/Platelet-interacting CAFs"                = "Erythroid+Platelet interacting CAFs"
##     )
##   - run_directional_ora_top10() and plot_dir_ora_top10() are defined above

## 7.1 Remap CAF cluster labels in DEG table
DEGs_Spatial_caf_dir <- DEGs_Spatial_caf_filtered %>%
  mutate(
    cluster_orig = cluster,
    cluster      = caf_cluster_remap[cluster_orig]
  ) %>%
  filter(!is.na(cluster))   # keep only CAF clusters that are in the remap

# Quick sanity check: how many genes per renamed CAF cluster?
DEGs_Spatial_caf_dir %>%
  count(cluster, name = "n_genes") %>%
  arrange(cluster)

## (Optional) Save the remapped DEG table for record-keeping
write_csv(
  DEGs_Spatial_caf_dir,
  "caf_deg_outputs/CAFs_DEGs_filtered_renamed_for_directional_ORA.csv"
)

## 7.2 Set output directory for Fig4ab CAF enrichment
out_dir_caf <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/CAF_enrichment_dir_top10"
dir.create(out_dir_caf, recursive = TRUE, showWarnings = FALSE)

## 7.3 Run directional ORA (Up vs Down) for CAF clusters
run_directional_ora_top10(
  deg_tbl = DEGs_Spatial_caf_dir,
  out_dir = out_dir_caf,
  prefix  = "CAF"
)

## 7.4 Generate KEGG + Hallmark top-10 barplots per CAF cluster
plot_dir_ora_top10(out_dir_caf, "CAF")





###############################################
## TAM ORA (KEGG + Hallmark) + Directional ORA
## With renamed TAM cluster labels
###############################################

## 1) Output paths for TAMs (Fig4ab)
tam_out_base <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/TAM_ORA"
tam_dir_out  <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/TAM_enrichment_dir_top10"

dir.create(tam_out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(tam_dir_out,  recursive = TRUE, showWarnings = FALSE)

## 2) Load TAM DEGs
##    Use the all-vs-rest filtered table (has cluster, gene, avg_log2FC)
if (!exists("DEGs_Spatial_tam_filtered")) {
  DEGs_Spatial_tam_filtered <- readr::read_csv(
    "tam_deg_outputs/TAM_DEGs__all_vs_rest__filtered.csv",
    show_col_types = FALSE
  )
}

stopifnot(all(c("cluster", "gene", "avg_log2FC") %in% colnames(DEGs_Spatial_tam_filtered)))

## 3) TAM cluster renaming (for ORA)
##    You can tweak these labels if you prefer slightly different wording.
tam_cluster_remap <- c(
  "M2-like/ECM-remodeling TAMs" = "IFN-TAM",
  "M1-like/Inflammatory TAMs"   = "proinflammatory TAM",
  "tTAMs"                       = "monocyte"
)

tam_levels <- c(
  "IFN-TAM",
  "proinflammatory TAM",
  "monocyte"
)

## keep original label, then overwrite `cluster` with renamed label
DEGs_Spatial_tam_filtered$cluster_original <- DEGs_Spatial_tam_filtered$cluster
DEGs_Spatial_tam_filtered$cluster <- tam_cluster_remap[DEGs_Spatial_tam_filtered$cluster]
DEGs_Spatial_tam_filtered$cluster <- factor(DEGs_Spatial_tam_filtered$cluster,
                                            levels = tam_levels)

# Quick sanity check
table(DEGs_Spatial_tam_filtered$cluster, useNA = "ifany")

###############################################
## 4) Standard ORA for TAMs (renamed clusters)
##    Uses the same run_enrichment_for_table() defined above for CAFs
###############################################

run_enrichment_for_table(
  deg_tbl = DEGs_Spatial_tam_filtered,
  out_dir = tam_out_base,
  prefix  = "TAM"
)

###############################################
## 5) Directional ORA (Up vs Down) for TAMs
##    (KEGG + Hallmark, using renamed TAM clusters)
###############################################

## 5.1 Create a TAM-specific DEG table for directional ORA
DEGs_Spatial_tam_dir <- DEGs_Spatial_tam_filtered %>%
  dplyr::mutate(
    cluster_orig = cluster,   # keep factor label (renamed)
    cluster      = as.character(cluster)  # make sure cluster is plain character
  ) %>%
  dplyr::filter(!is.na(cluster))

# Optional sanity check: number of genes per renamed TAM cluster
DEGs_Spatial_tam_dir %>%
  dplyr::count(cluster, name = "n_genes") %>%
  dplyr::arrange(cluster)

# Optional: save this table
readr::write_csv(
  DEGs_Spatial_tam_dir,
  "tam_deg_outputs/TAM_DEGs__all_vs_rest__filtered_renamed_for_directional_ORA.csv"
)

## 5.2 Run directional ORA for TAMs
run_directional_ora_top10(
  deg_tbl = DEGs_Spatial_tam_dir,
  out_dir = tam_dir_out,
  prefix  = "TAM"
)

## 5.3 Generate KEGG + Hallmark top-10 barplots per TAM cluster
plot_dir_ora_top10(tam_dir_out, "TAM")




###############################################
## B cell ORA (KEGG + Hallmark) + Directional ORA
## With renamed B cell cluster labels
###############################################


## 1) Output paths for B cells (Fig4ab)
bcell_out_base <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/Bcell_ORA"
bcell_dir_out  <- "<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Visualization/Fig4ab/Bcell_enrichment_dir_top10"

dir.create(bcell_out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(bcell_dir_out,  recursive = TRUE, showWarnings = FALSE)

## 2) Load B cell DEGs
##    ⚠️ Adjust the file path below if your B cell DEG table
##    has a different name or lives in another folder.
if (!exists("DEGs_Spatial_bcell_filtered")) {
  DEGs_Spatial_bcell_filtered <- readr::read_csv(
    "bcell_deg_outputs/Bcell_DEGs__all_vs_rest__filtered.csv",
    show_col_types = FALSE
  )
}

stopifnot(all(c("cluster", "gene", "avg_log2FC") %in% colnames(DEGs_Spatial_bcell_filtered)))

## 3) B cell cluster renaming (for ORA)
##    You can tweak these labels if you prefer slightly different text.
bcell_cluster_remap <- c(
  "tBcells"                   = "tB cells",
  "Mature B cell/Plasma cell" = "Mature/Plasma B cells",
  "Mucosal Plasma B cells"    = "Mucosal plasma B cells"
)

bcell_levels <- c(
  "tB cells",
  "Mature/Plasma B cells",
  "Mucosal plasma B cells"
)

DEGs_Spatial_bcell_filtered$cluster_original <- DEGs_Spatial_bcell_filtered$cluster
DEGs_Spatial_bcell_filtered$cluster <- bcell_cluster_remap[DEGs_Spatial_bcell_filtered$cluster]
DEGs_Spatial_bcell_filtered$cluster <- factor(DEGs_Spatial_bcell_filtered$cluster,
                                              levels = bcell_levels)

# Sanity check
table(DEGs_Spatial_bcell_filtered$cluster, useNA = "ifany")

###############################################
## 4) Standard ORA for B cells (renamed clusters)
###############################################

run_enrichment_for_table(
  deg_tbl = DEGs_Spatial_bcell_filtered,
  out_dir = bcell_out_base,
  prefix  = "Bcell"
)

###############################################
## 5) Directional ORA (Up vs Down) for B cells
###############################################

## 5.1 Create B cell DEG table for directional ORA
DEGs_Spatial_bcell_dir <- DEGs_Spatial_bcell_filtered %>%
  dplyr::mutate(
    cluster_orig = cluster,
    cluster      = as.character(cluster)  # ensure simple character labels
  ) %>%
  dplyr::filter(!is.na(cluster))

# Optional: sanity check
DEGs_Spatial_bcell_dir %>%
  dplyr::count(cluster, name = "n_genes") %>%
  dplyr::arrange(cluster)

# Optional: save remapped table
readr::write_csv(
  DEGs_Spatial_bcell_dir,
  "bcell_deg_outputs/Bcell_DEGs__all_vs_rest__filtered_renamed_for_directional_ORA.csv"
)

## 5.2 Run directional ORA for B cells
run_directional_ora_top10(
  deg_tbl = DEGs_Spatial_bcell_dir,
  out_dir = bcell_dir_out,
  prefix  = "Bcell"
)

## 5.3 Generate KEGG + Hallmark top-10 barplots per B cell cluster
plot_dir_ora_top10(bcell_dir_out, "Bcell")
