# Part 4 — Visualization (publication-ready)

This folder contains **figure-only** scripts used to generate the manuscript figures.
Heavy computation (QC, integration, clustering, inferCNV/CellChat running) should be done in Parts 1–3.
Part 4 assumes those outputs exist and **only reads processed objects** (RDS/CSV) to generate plots.

## Files

- `00_common.R`  
  Shared helpers used by all figure scripts:
  - package loading (minimal set)
  - a single `theme_pub()` applied globally
  - `save_plot()` wrapper that standardizes PDF/PNG export

- `fig1.R` — Cohort overview (spatial + UMAP)
- `fig2.R` — Malignant signatures / subclustering panels (no biological renaming)
- `fig3.R` — Additional malignant characterization panels
- `fig5.R` — Additional analysis panels (composition, fibro/macro panels, etc.)
- `suppl_fig1.R`, `suppl_fig3.R`, `suppl_fig4.R` — Supplementary figures

> **Note:** `fig4.R` was not included in this refinement bundle (not uploaded in this chat). Re-upload it to refine it the same way.

## How to run

From the directory containing the figure scripts:

```bash
# (recommended) start a clean R session
R

# in R:
source("fig1.R")
source("fig2.R")
source("fig3.R")
source("fig5.R")
source("suppl_fig1.R")
source("suppl_fig3.R")
source("suppl_fig4.R")
```

Each script writes outputs under `figures/<figure_name>/` by default.

## What you must edit

At the top of each script there is an **Inputs** block (paths like `outputs/part2/...`).
Replace these with the correct paths on your machine or in your project structure.

Any remaining placeholder strings like `<LOCAL_PATH>/...` are intentional markers that you must replace.

## Style conventions

- No `setwd()` anywhere (use `out_dir` + `file.path()`).
- No `install.packages()` / `BiocManager::install()` inside scripts.
- One global theme (`theme_set(theme_pub())`) for consistent publication styling.
- Use `save_plot()` / `ggsave()` with explicit width/height.
