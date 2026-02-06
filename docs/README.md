# Spatial Analysis HNSCC Framework

**Authors:** [Gregoire Marret](https://github.com/gmarret), [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA), [Jinsu An](https://github.com/Jinsuan-UofC)


**Contact:** [gregoire.marret@uhn.ca](gregoire.marret@uhn.ca), [farnoosh.abbasaghababazadeh@uhn.ca](farnoosh.abbasaghababazadeh@uhn.ca), [jinsu.an@ucalgary.ca](jinsu.an@ucalgary.ca)

**Description:** A Spatial Analysis Framework for Head and Neck Squamous Cell Carcinoma

This repository contains computational workflows for analyzing spatial transcriptomics data from recurrent/metastatic head and neck squamous cell carcinoma (RM-HNSCC). The analyses characterize spatial tumor–immune–stromal organization, infer cell–cell communication, and derive spatially informed gene expression signatures that predict response to immune checkpoint blockade across cancer types.

--------------------------------------

[![pixi-badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json&style=flat-square)](https://github.com/prefix-dev/pixi)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json&style=flat-square)](https://github.com/astral-sh/ruff)
[![Built with Material for MkDocs](https://img.shields.io/badge/mkdocs--material-gray?logo=materialformkdocs&style=flat-square)](https://github.com/squidfunk/mkdocs-material)

![GitHub last commit](https://img.shields.io/github/last-commit/bhklab/spatial-analysis-framework-hnscc?style=flat-square)
![GitHub issues](https://img.shields.io/github/issues/bhklab/spatial-analysis-framework-hnscc?style=flat-square)
![GitHub pull requests](https://img.shields.io/github/issues-pr/bhklab/spatial-analysis-framework-hnscc?style=flat-square)
![GitHub contributors](https://img.shields.io/github/contributors/bhklab/spatial-analysis-framework-hnscc?style=flat-square)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/bhklab/spatial-analysis-framework-hnscc?style=flat-square)

## Project Overview

This repository serves as the computational companion to a spatial transcriptomics study of RM-HNSCC treated with immune checkpoint blockade. It implements the analytical workflow underlying the manuscript, organizing the full sequence of computational steps required to reproduce the reported results given appropriate data access.

The framework supports:

- preprocessing, integration, and quality control of spatial transcriptomics data;
- annotation of malignant and non-malignant compartments using orthogonal approaches;
- modeling of spatial neighborhoods and intercellular communication;
- extraction of spatial tumor–microenvironment features; and
- translation of spatial features into gene expression signatures validated in external bulk RNA-seq cohorts.

All analyses are implemented as modular, configurable scripts designed for transparency, reproducibility, and extensibility. Protected patient-level data are not distributed in this repository; however, all code required to reproduce the analyses is provided.

**Reproducibility**

- Modular pipeline implemented in **R** and **Python**
- Dependency management via **Pixi**
- Configuration-driven execution to support multiple cohorts and analyses

---

## Set Up

### Prerequisites

Pixi is required to run this project.
If you haven't installed it yet, [follow these instructions](https://pixi.sh/latest/)

### Installation

```bash
# Clone the repository
git clone https://github.com/bhklab/spatial-analysis-framework-hnscc.git
cd spatial-analysis-framework-hnscc

# Install dependencies via Pixi
pixi install
```

## Repository Structure

```
spatial-analysis-framework-hnscc/
├── config/           # YAML config files for each dataset and center
├── data/             # Raw, processed, and results directories
├── workflow/         # Scripts for analysis
├── docs/             # MkDocs-based project documentation
│   └── README.md     # Documentation index and setup instructions
└── pixi.toml         # Pixi environment specification
```

---

## Documentation

Full documentation, including usage instructions, data setup, config templates, and pipeline stages, will be available in the `docs/` folder or via published GitHub Pages.

Start by downloading and organizing the raw input datasets as described in [`data/rawdata/README.md`](https://github.com/bhklab/spatial-analysis-framework-hnscc/blob/main/data/rawdata/README.md).

For data download and processing, please refer to the univariable repository:  
🔗 [https://github.com/bhklab/spatial-analysis-framework-hnscc?tab=readme-ov-file](https://github.com/bhklab/spatial-analysis-framework-hnscc?tab=readme-ov-file)
