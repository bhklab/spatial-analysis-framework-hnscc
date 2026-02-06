# -----------------------------------------------------------------------------------------------------
# Description:
# This script performs meta-analyses of gene signature associations with immunotherapy outcomes
# across multiple independent cohorts. It integrates study-level association results derived
# from precomputed regression models and summarizes effects using fixed- or random-effects
# meta-analysis via the `metafun()` function.
#
# The meta-analyses are conducted at multiple levels:
#   - Pan-cancer (all eligible cancer types combined)
#   - Cancer-specific (e.g., Melanoma, Lung, Kidney)
#   - Treatment-specific (e.g., PD-(L)1)
#
# Analyses include:
#   - Progression-free survival (PFS) associations from Cox proportional hazards models
#     using either continuous signature metrics (e.g., distance to centroid2) or
#     dichotomized high/low groupings
#   - Objective response (responder vs. non-responder) associations from logistic
#     regression models using continuous signature metrics
#
# Input:
#   - Precomputed association result files in RDA format located in `data/results/assoc/`:
#       * sig_pfs.rda        — Cox models (continuous predictors)
#       * sig_pfs_dicho.rda  — Cox models (dichotomized predictors)
#       * sig_logreg.rda     — Logistic regression models (response)
#   - Each file must contain study-level estimates including:
#       Gene, Coef, SE, Pval, Study, N, Cancer_type, Treatment
#
# Output:
#   - Meta-analysis summary tables written to `data/results/meta/`, including:
#       * Pooled effect sizes and confidence intervals
#       * P-values and heterogeneity metrics (I², Q-test)
#       * Stratification by outcome (PFS, R/NR), predictor type (continuous, dichotomized),
#         and analysis group (pan-cancer or cancer-specific)
#
# Notes:
#   - Only features with results from ≥3 independent studies are included
#   - Meta-analyses are performed using study-level effect sizes and standard errors
# -----------------------------------------------------------------------------------------------------
##################################################
## Load libraries
##################################################

library(meta)
library(forestplot)
library(dplyr)
library(data.table)
library(PredictioR)
library(MultiAssayExperiment)

##################################################
## Set up directory
##################################################

dir.input <- 'data/results/assoc'
dir.output <- 'data/results/meta'

########################################################################################################
############################################## pan-cancer (PFS) ########################################
########################################################################################################
## Continuous signature ----> PFS
load(file.path(dir.input, "sig_pfs.rda"))
genes <- unique(res_pfs$Gene)

res_meta <- lapply(1:length(genes), function(k){

  df <- res_pfs[res_pfs$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Cont"
res_meta$group <- "Pan"
res_meta$Outcome <- "PFS"

write.csv(res_meta, file = file.path(dir.output, 'meta_pancancer_pfs_cont.csv'), row.names= FALSE)

## Dicho signature ---> PFS
load(file.path(dir.input, "sig_pfs_dicho.rda"))
genes <- unique(res_pfs$Gene)

res_meta <- lapply(1:length(genes), function(k){

  df <- res_pfs[res_pfs$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Dicho"
res_meta$group <- "Pan"
res_meta$Outcome <- "PFS"

write.csv(res_meta, file = file.path(dir.output, 'meta_pancancer_pfs_dicho.csv'), row.names= FALSE)

## Continuous ---> R/NR 
load(file.path(dir.input, "sig_logreg.rda"))
genes <- unique(res_logreg$Gene)

res_meta <- lapply(1:length(genes), function(k){

  df <- res_logreg[res_logreg$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Cont"
res_meta$group <- "Pan"
res_meta$Outcome <- "R/NR"

write.csv(res_meta, file = file.path(dir.output, 'meta_pancancer_logreg.csv'), row.names= FALSE)

########################################################################################################
############################################## per-cancer (Melanoma) ###################################
########################################################################################################
## Continuous ---> R/NR  

load(file.path(dir.input, "sig_logreg.rda"))
res_logreg <- res_logreg[res_logreg$Cancer_type == 'Melanoma', ]
genes <- unique(res_logreg$Gene)

res_meta <- lapply(1:length(genes), function(k){
  
  df <- res_logreg[res_logreg$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Cont"
res_meta$group <- "Melanoma"
res_meta$Outcome <- "R/NR"

write.csv(res_meta, file = file.path(dir.output, 'meta_melanoma_logreg.csv'), row.names= FALSE)

########################################################################################################
############################################## per-cancer (Lung) #######################################
########################################################################################################
## Continuous signature ----> PFS 
load(file.path(dir.input, "sig_pfs.rda"))
res_pfs <- res_pfs[res_pfs$Cancer_type == 'Lung', ]
genes <- unique(res_pfs$Gene)

res_meta <- lapply(1:length(genes), function(k){
  
  df <- res_pfs[res_pfs$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta  <- res_meta [!is.na(res_meta $Coef), ]
res_meta $Type <- "Cont"
res_meta $group <- "Lung"
res_meta$Outcome <- "PFS"

write.csv(res_meta, file = file.path(dir.output, 'meta_lung_pfs_cont.csv'), row.names= FALSE)

## Dicho signature ----> PFS
load(file.path(dir.input, "sig_pfs_dicho.rda"))
res_pfs <- res_pfs[res_pfs$Cancer_type == 'Lung', ]
genes <- unique(res_pfs$Gene)

res_meta <- lapply(1:length(genes), function(k){

  
  df <- res_pfs[res_pfs$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta $Coef), ]
res_meta  $Type <- "Dicho"
res_meta $group <- "Lung"
res_meta$Outcome <- "PFS"

write.csv(res_meta, file = file.path(dir.output, 'meta_lung_pfs_dicho.csv'), row.names= FALSE)

## Continuous ---> R/NR  

load(file.path(dir.input, "sig_logreg.rda"))
res_logreg <- res_logreg[res_logreg$Cancer_type == 'Lung', ]
genes <- unique(res_logreg$Gene)

res_meta <- lapply(1:length(genes), function(k){

  df <- res_logreg[res_logreg$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Cont"
res_meta$group <- "Lung"
res_meta$Outcome <- "R/NR"

write.csv(res_meta, file = file.path(dir.output, 'meta_lung_logreg.csv'), row.names= FALSE)

########################################################################################################
############################################## per-cancer (Kidney) #####################################
########################################################################################################
## Continuous ---> R/NR  
load(file.path(dir.input, "sig_logreg.rda"))
res_logreg <- res_logreg[res_logreg$Cancer_type == 'Kidney', ]
genes <- unique(res_logreg$Gene)

res_meta <- lapply(1:length(genes), function(k){

  df <- res_logreg[res_logreg$Gene == genes[k], ]
  if(nrow(df) >= 3){

   res <- metafun(coef = df$Coef,
                     se = df$SE,
                     study  = df$Study,
                     pval = df$Pval,
                     n = df$N,
                     cancer.type = df$Cancer_type,
                     treatment = df$Treatment,
                     feature = unique(df$Gene),
                     cancer.spec = FALSE,
                     treatment.spec = FALSE)
   res$meta_summery

  }else{

   data.frame(Gene = genes[k],
              Coef = NA,
              SE =NA,
              CI_lower = NA,
              CI_upper = NA,
              Pval = NA,
              I2 = NA,
              Q_Pval = NA)

  }

})

res_meta <- do.call(rbind, res_meta)        
res_meta <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Cont"
res_meta$group <- "Kidney"
res_meta$Outcome <- "R/NR"

write.csv(res_meta, file = file.path(dir.output, 'meta_kidney_logreg.csv'), row.names= FALSE)

###################################################
## Merge integration results
###################################################
meta.files <- list.files(dir.output)
meta_res <- lapply(1:length(meta.files), function(k){

  read.csv(file.path(dir.output, meta.files[k]))

})

meta_res <- do.call(rbind, meta_res)
write.csv(meta_res, file = file.path(dir.output, 'meta_all.csv'), row.names = FALSE)
