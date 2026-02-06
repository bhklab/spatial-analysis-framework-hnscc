# -----------------------------------------------------------
# Description:
# This script reads input and output directory paths and study metadata
# from `config_proc.yaml`, and performs the following:
#
#   - Loads raw MultiAssayExperiment (MAE) object from `data/rawdata`
#   - Converts MAE to SummarizedExperiment (SE)
#   - Loads gene signature data and metadata
#
# Config file: config/config_proc.yaml
# -----------------------------------------------------------
########################################################
## Load Libraries
########################################################
library(MultiAssayExperiment)
library(data.table)
library(PredictioR)
library(yaml)

###########################################################
## Set up working directory and study configuration
###########################################################
# Load configuration file
config <- yaml::read_yaml("config/config_proc.yaml")

dir_in <- config$dir_in # "data/rawdata"
dir_out <- config$dir_out # "data/procdata"

study_icb <- config$study_icb # "ICB_Gide"
cancer_type <- config$cancer_type # "Melanoma"
treatment_type <- config$treatment_type # "PD-(L)1"  (Other options include: CTLA-4, IO+combo, etc.)

########################################################
## Load Raw ICB Dataset
########################################################
# Ensure the following file is downloaded from:
# https://www.orcestra.ca/clinical_icb/62f29e85be1b2e72a9c177f4
# and placed in `data/rawdata/ICB_Gide.rds`

mae_obj <- readRDS(file.path(dir_in, "ICB_Gide.rds"))  # MultiAssayExperiment
se_obj  <- PredictioR::createSE(mae_obj)  # Extract SummarizedExperiment (TPM + clinical)

