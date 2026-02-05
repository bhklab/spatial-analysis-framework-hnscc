# -----------------------------------------------------------------------------------------------------
# Description:
# This script performs meta-analysis of gene signature associations with immunotherapy outcomes 
# across multiple studies. It combines results from precomputed association tests (e.g., Cox or 
# logistic regression) and applies fixed/random effects meta-analysis using the `metafun()` function.
#
# The analyses cover:
#   - Pan-cancer meta-analysis (across all cancer types)
#   - Per-cancer meta-analysis (e.g., Melanoma, Lung)
#   - Per-treatment meta-analysis (e.g., PD-(L)1)
#
# For each category, both continuous and dichotomized (high/low) signature scores are analyzed:
#   - PFS outcomes (Cox regression)
#   - Response vs. Non-response outcomes (Logistic regression)
#
# Input:
#   - Precomputed signature association results in RDA format:
#       * sig_pfs.rda: Cox regression on continuous scores
#       * sig_pfs_dicho.rda: Cox regression on dichotomized scores
#       * sig_logreg.rda: Logistic regression on response (R/NR)
#   - Each file must contain variables: Gene, Coef, SE, Pval, Study, N, Cancer_type, Treatment
#
# Output:
#   - Meta-analysis summary tables saved as CSV files in: result/meta/
#       * Includes effect size, confidence interval, p-value, heterogeneity metrics (I², Q-p)
#       * Results categorized by analysis type (Cont/Dicho), outcome type (PFS/Logreg), and group
#
# Dependencies:
#   - Libraries: meta, forestplot, dplyr, data.table, PredictioR, MultiAssayExperiment
#   - Custom function: metafun() (assumed to be defined in the environment)
#
# Notes:
#   - Only genes with results from ≥3 studies are included in the meta-analysis.
#   - Meta-analysis is performed using study-level effect sizes and standard errors.
#
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

dir.input <- 'result/assoc'
dir.output <- 'result/meta'

########################################################################################################
############################################## pan-cancer (PFS) ########################################
########################################################################################################
## Continuous signature ----> PFS
load(file.path(dir.input, "sig_pfs.rda"))
genes <- unique(res_pfs$Gene)
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
## Continuous signature ---> PFS
load(file.path(dir.input, "sig_pfs.rda"))
res_pfs <- res_pfs[res_pfs$Cancer_type == 'Melanoma', ]
genes <- unique(res_pfs$Gene)
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
res_meta  <- res_meta[!is.na(res_meta$Coef), ]
res_meta$Type <- "Cont"
res_meta $group <- "Melanoma"
res_meta$Outcome <- "PFS"

write.csv(res_meta, file = file.path(dir.output, 'meta_melanoma_pfs_cont.csv'), row.names= FALSE)

## Dicho signature ---> PFS
load(file.path(dir.input, "sig_pfs_dicho.rda"))
res_pfs <- res_pfs[res_pfs$Cancer_type == 'Melanoma', ]

genes <- unique(res_pfs$Gene)
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
res_meta  <- res_meta [!is.na(res_meta$Coef), ]
res_meta$Type <- "Dicho"
res_meta$group <- "Melanoma"
res_meta$Outcome <- "PFS"

write.csv(res_meta, file = file.path(dir.output, 'meta_melanoma_pfs_dicho.csv'), row.names= FALSE)

## Continuous ---> R/NR  

load(file.path(dir.input, "sig_logreg.rda"))
res_logreg <- res_logreg[res_logreg$Cancer_type == 'Melanoma', ]
genes <- unique(res_logreg$Gene)
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]
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
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
genes <- genes[!genes %in% c('CIN25_Carter', 'CIN70_Carter')]

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
