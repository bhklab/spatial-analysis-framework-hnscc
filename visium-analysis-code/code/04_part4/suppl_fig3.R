# ============================================================
# Supplementary Figure 3
# Publication polished script (paths parameterized, no setwd)
# ============================================================


# ---- shared helpers ----
# Assumes this script lives in the same folder as 00_common.R.
# If not, change the path below (e.g., source("code/04_part4/00_common.R")).
source("00_common.R")
theme_set(theme_pub())

suppressPackageStartupMessages({
  library(forcats)
  library(ggplot2)
  library(stringr)
})

# ----------------------------
# Inputs (edit these paths)
# ----------------------------


# ----------------------------
# Outputs
# ----------------------------
out_dir <- file.path("figures", "suppl_fig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



#IPA 
## ==========================================================
## 0. Setup
## ==========================================================
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part2/2.3 - Visualization/IPA/")  # removed (publication): use out_dir + file.path()

## All IPA txt exports (named by subcluster)
ipa_files <- c(
  'cycling' = "ipa_Cancer_C0.txt",
  'tumor core' = "ipa_Cancer_C1.txt",
  'neutrophil-inflamed' = "ipa_Cancer_C2.txt",
  'fibrovascular niche' = "ipa_Cancer_C3.txt",
  'p-EMT' = "ipa_Cancer_C4.txt"
)

## ==========================================================
## 1. Import IPA canonical table
## ==========================================================
importIPAenrichment <- function(file) {
  df <- read.delim(
    file,
    header = TRUE,
    skip = 24,
    check.names = FALSE,
    fill = TRUE,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  df <- df[, 1:5]
  colnames(df) <- c("Pathway", "-log(p-value)", "zScore", "Ratio", "Molecules")
  
  # Clean / trim pathway strings
  df$Pathway <- trimws(as.character(df$Pathway))
  
  # Remove junk rows, including literal "NA"/"N/A" etc and numeric-only rows
  bad_na_vals <- c("NA", "Na", "na", "N/A", "n/a")
  df <- df[
    !is.na(df$Pathway) &
      df$Pathway != "" &
      !(df$Pathway %in% bad_na_vals) &
      !grepl("^[0-9.]+$", df$Pathway),
  ]
  
  df
}



## Load all IPA files into a named list
ipa_l <- lapply(ipa_files, importIPAenrichment)
topPathways <- lapply(ipa_l, get_top_pathways)
names(topPathways) <- names(ipa_files)

cluster_names <- gsub("Cancer_C", "Cluster ", names(topPathways))
names(topPathways) <- cluster_names


## ==========================================================
## 2. Get top activated & inhibited pathways per file
## ==========================================================
get_top_pathways <- function(df, n = 10, p_cutoff = 0.05) {
  df$zScore <- suppressWarnings(as.numeric(df$zScore))
  df$logp   <- suppressWarnings(as.numeric(df$`-log(p-value)`))
  df$pvalue <- 10^(-df$logp)
  
  # remove empty names again just to be safe
  df <- df[!is.na(df$Pathway) & trimws(df$Pathway) != "", ]
  
  sig <- subset(df, !is.na(pvalue) & pvalue < p_cutoff & !is.na(zScore))
  
  top_activated <- head(sig[order(-sig$zScore), ], n)
  top_inhibited <- head(sig[order(sig$zScore), ], n)
  
  list(
    top_activated = top_activated,
    top_inhibited = top_inhibited
  )
}

topPathways <- lapply(ipa_l, get_top_pathways)
names(topPathways) <- names(ipa_files)
plot_top_pathways_bar <- function(top_pathways, filename_prefix) {
  # Combine activated & inhibited, mark direction
  plotdf <- rbind(
    within(top_pathways$top_activated, { Direction <- "Activated" }),
    within(top_pathways$top_inhibited, { Direction <- "Inhibited" })
  )
  
  # Clean / trim pathway names
  plotdf$Pathway <- trimws(as.character(plotdf$Pathway))
  
  # Drop obvious junk
  bad_na_vals <- c("NA", "Na", "na", "N/A", "n/a")
  plotdf <- plotdf[!duplicated(plotdf$Pathway), ]
  plotdf <- subset(
    plotdf,
    !is.na(Pathway) &
      Pathway != "" &
      !(Pathway %in% bad_na_vals) &
      !is.na(Direction)
  )
  
  # Original and absolute z for plotting
  plotdf$zScore_orig <- suppressWarnings(as.numeric(plotdf$zScore))
  plotdf$zScore_plot <- abs(plotdf$zScore_orig)
  
  # Wrap pathway names
  plotdf$Pathway_wrapped <- stringr::str_wrap(plotdf$Pathway, width = 35)
  
  # *** CRITICAL CLEANING STEP ***
  # Drop any rows where the wrapped label or z is NA or literal "NA"
  plotdf <- subset(
    plotdf,
    !is.na(Pathway_wrapped) &
      as.character(Pathway_wrapped) != "NA" &
      !is.na(zScore_plot)
  )
  
  # Order: activated (high z) first, inhibited (most negative) next
  plotdf <- rbind(
    plotdf[plotdf$Direction == "Activated", ][order(-plotdf$zScore_orig), ],
    plotdf[plotdf$Direction == "Inhibited", ][order(plotdf$zScore_orig), ]
  )
  
  # Final factor levels and droplevels to strip any hidden NA level
  plotdf$Pathway_wrapped <- droplevels(
    factor(plotdf$Pathway_wrapped,
           levels = rev(unique(plotdf$Pathway_wrapped)))
  )
  
  # p-value labels
  plotdf$pval_lab <- ifelse(
    is.na(plotdf$pvalue), "",
    paste0("p=", formatC(plotdf$pvalue, format = "e", digits = 2))
  )
  
  max_z <- max(plotdf$zScore_plot, na.rm = TRUE)
  
  p <- ggplot(plotdf, aes(x = zScore_plot, y = Pathway_wrapped, fill = Direction)) +
    geom_col(width = 0.75, color = "black", alpha = 0.9) +
    geom_text(aes(label = pval_lab), hjust = -0.1, size = 5) +
    scale_fill_manual(
      values = c("Activated" = "#D73027", "Inhibited" = "#4575B4"),
      name = NULL,
      na.translate = FALSE
    ) +
    scale_y_discrete(drop = TRUE, na.translate = FALSE) +   # <- tell ggplot to drop NA on y
    labs(
      x = "absolute IPA z-score",
      y = NULL,
      subtitle = filename_prefix
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.2)),
      limits = c(0, max_z * 1.2)
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "top",
      axis.text.y     = element_text(size = 14, face = "bold"),
      axis.text.x     = element_text(size = 12),
      axis.title.x    = element_text(size = 14, face = "bold"),
      plot.subtitle   = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.ticks.y    = element_blank(),
      legend.text     = element_text(size = 13),
      plot.margin     = margin(12, 12, 12, 12)
    )
  
  ggsave(
    paste0(filename_prefix, "_top_pathways_barplot_absZ.png"),
    plot = p,
    width = 8,
    height = max(6, nrow(plotdf) * 0.35),
    dpi = 600
  )
  
  p
}
# setwd("<LOCAL_PATH>/jinsu/Desktop/UHN - Lillian Siu/INSPIRE and IO-kin/v3/Part4_Revision_v2/Supp_Fig3/")  # removed (publication): use out_dir + file.path()
for (nm in names(topPathways)) {
  plot_top_pathways_bar(topPathways[[nm]], nm)
}