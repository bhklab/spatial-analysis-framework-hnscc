# ============================================================
# 00_common.R — shared helpers for Part 4 (Visualization)
# - Centralizes packages, plotting theme, and save helpers
# - Keeps figure scripts short and non-repetitive
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
})

# ---- global ggplot theme (edit once; applied everywhere) ----
theme_pub <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "plain"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

# ---- I/O helpers ----
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

assert_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  invisible(path)
}

# A single save wrapper used across scripts
save_plot <- function(p, filename, width, height, dpi = 600, ...) {
  ext <- tolower(tools::file_ext(filename))
  ensure_dir(dirname(filename))

  if (ext == "pdf") {
    ggsave(filename = filename, plot = p, width = width, height = height,
           units = "in", device = cairo_pdf, ...)
  } else if (ext %in% c("png", "tiff", "jpeg", "jpg")) {
    ggsave(filename = filename, plot = p, width = width, height = height,
           units = "in", dpi = dpi, ...)
  } else {
    # fall back to default
    ggsave(filename = filename, plot = p, width = width, height = height,
           units = "in", dpi = dpi, ...)
  }
  invisible(filename)
}
