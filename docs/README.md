# Spatial Analysis HNSCC Framework

**Authors:** [Gregoire Marret](https://github.com/gmarret), [Jinsu An](https://github.com/Jinsuan-UofC), [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA)


**Contact:** [gregoire.marret@uhn.ca](gregoire.marret@uhn.ca), [jinsu.an@ucalgary.ca](jinsu.an@ucalgary.ca), [farnoosh.abbasaghababazadeh@uhn.ca](farnoosh.abbasaghababazadeh@uhn.ca)

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

This repository implements ....

- **Reproducibility**: modular pipeline using **Pixi**, **Python**, and **R**

---

## Set Up

### Prerequisites

Pixi is required to run this project.
If you haven't installed it yet, [follow these instructions](https://pixi.sh/latest/)

### Installation

```bash
# Clone the repository
git clone https://github.com/bhklab/predictio-uv-dist.git
cd predictio-uv-dist

# Install dependencies via Pixi
pixi install
```

## Repository Structure

```
predictio-uv-dist/
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
