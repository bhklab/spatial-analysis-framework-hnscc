# -------------------------------------------------------------------------------
# Description:
# This script generates visualizations to assess the association between a
# spatially informed gene expression signature (SpatialSignature_Marret),
# quantified as a distance-based metric, and immunotherapy outcomes across
# multiple bulk RNA-seq cohorts.
#
# The signature represents the distance of each sample to a reference centroid
# associated with favorable immunotherapy outcomes and is evaluated against
# progression-free survival (PFS) and objective response (R vs. NR).
#
# Visualizations generated include:
#
# 1. Kaplan–Meier survival curves for PFS using median-based dichotomization
#    of the spatial signature within each cohort.
#
# 2. Boxplots comparing continuous spatial signature values between responders
#    and non-responders, with Wilcoxon rank-sum test p-values.
#
# 3. Forest plots summarizing study-level association results derived from:
#    - Cox proportional hazards models for PFS (logHR)
#    - Logistic regression models for response (logOR)
#    Stratified by:
#      • Pan-cancer
#      • Cancer type (e.g., Lung, Melanoma, Kidney)
#      • Treatment category (e.g., PD-(L)1)
#
# 4. Meta-analysis–based visual summaries using pooled effect estimates
#    obtained from external association analyses.
#
# Required Inputs:
#   - Precomputed spatial signature values from `data/results/score/`
#   - Clinical outcome data from `data/results/data/`
#   - Association result tables from `data/results/assoc/`
#
# Output:
#   - Kaplan–Meier plots saved to `data/results/Fig/KM_median/`
#   - Boxplots saved to `data/results/Fig/boxplot/`
#   - Forest plots saved to `data/results/Fig/forestplot/`
#
# Notes:
#   - PFS is administratively censored at 24 months.
#   - Median-based dichotomization requires ≥3 samples per group.
#   - Forest plots display log hazard ratios (PFS) or log odds ratios (response)
#     with corresponding confidence intervals.
# -------------------------------------------------------------------------------
############################################################################
# Load libraries
############################################################################
library(stringr)
library(meta)
library(survival)
library(survminer)
library(forestplot)
library(dplyr)
library(data.table)
library(PredictioR)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(MultiAssayExperiment)

dir.score <- 'data/results/score'
dir.input <- 'data/results/assoc'
dir.output <- 'data/results/Fig'

sig_name <- 'SpatialSignature_Marret'
####################################################################
## KM plot --> median cut-off
####################################################################
## load signatures
files.sig <- list.files(dir.score)

geneSig.score <- lapply(1:length(files.sig), function(k){
  load(file.path(dir.score, files.sig[k]))
  geneSig.score
})

study_icb <- substr(files.sig, 1, nchar(files.sig) - 4)
names(geneSig.score) <- study_icb

## load data
bin.cutoff <- 0.50
files <- list.files(file.path('data/results/data'))

for(k in 1:length(geneSig.score)){
  
  print(files[k])
  load(file.path('data/results/data', files[k]))

  var.pfs <- colData(dat_icb)$event_occurred_pfs

  if(sum(!is.na(var.pfs)) > 0 & sum(var.pfs != 0, na.rm =TRUE) > 0 ){

  df <- geneSig.score[[k]]
  df <- df[rownames(df) == sig_name, ]
  study_name <- names(geneSig.score)[k]
  study_icb <- strsplit(study_name, '__')[[1]][1]
  cancer_type <- strsplit(study_name, '__')[[1]][2]
  treatment_type <- strsplit(study_name, '__')[[1]][3]
  
pdf(file.path(dir.output, 'KM_median', paste(names(geneSig.score)[k], 'KM.pdf', sep='__')), height = 5, width = 5)

KMPlot(status = dat_icb$event_occurred_pfs,
       time = dat_icb$survival_time_pfs,
       time.censor = 24,
       var =  as.numeric(df),
       title = " ",
       xlab = "Time (Months)",
       ylab = "Progression Free Survival", 
       method = "median",
       n0.cutoff = 3, 
       n1.cutoff = 3,
       var.type = FALSE)

dev.off()

  }else{
     print('no PFS in a given cohort')
  }

}


####################################################################
## Boxplot plot
####################################################################
## load signatures
files.sig <- list.files(dir.score)

geneSig.score <- lapply(1:length(files.sig), function(k){
  load(file.path(dir.score, files.sig[k]))
  geneSig.score
})

study_icb <- substr(files.sig, 1, nchar(files.sig) - 4)
names(geneSig.score) <- study_icb

## load data
files <- list.files(file.path('data/results/data'))

for(k in 1:length(geneSig.score)){
  
  load(file.path('data/results/data', files[k]))
  
  print(files[k])
  clin <- colData(dat_icb)
  df <- geneSig.score[[k]]
  if(sum(rownames(df) == sig_name) != 0){
  
  df <- as.numeric(df[rownames(df) == sig_name, ])
  clin$sig <- df
  clin <- clin[!is.na(clin$response), ]
  
    if(sum(!is.na(clin$response)) > 0 & sum(clin$response == 'R', na.rm =TRUE) > 2 ){

  study_name <- names(geneSig.score)[k]
  study_icb <- strsplit(study_name, '__')[[1]][1]
  cancer_type <- strsplit(study_name, '__')[[1]][2]
  treatment_type <- strsplit(study_name, '__')[[1]][3]
  
pdf(file.path(dir.output, 'boxplot', paste(names(geneSig.score)[k], 'boxplot.pdf', sep='__')), 
    height = 4, width = 4)

p <- ggplot(clin, aes(x = response, y = sig, fill = response)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.format") + 
  scale_fill_manual(values = c("R" = "#1b9e77", "NR" = "#7570b3")) +
  labs(y = "signature score", x = "") +
   theme(
      axis.text.x=element_text(size=10,  face="bold"),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y=element_text(size=10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position="bottom",
      legend.text = element_text(size = 10, face="bold"),
      legend.title = element_blank())

print(p)

dev.off()

  }else{
     print('not min 3 response')
  }

  }else{

     print('signature not available')
  }
  
}

################################################
## Forestplot: PFS Pan
################################################

df <- read.csv(file = file.path(dir.input, "sig_pfs.csv"))
df <- df[df$Gene == sig_name, ]

pdf(file=file.path(dir.output, 'forestplot/PFS',
                    paste("INSPIRE_IOKIN", "pan_cont.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logHR estimate",
                 label = "logHR")

dev.off()


df <- read.csv(file = file.path(dir.input, "sig_pfs_dicho.csv"))
df <- df[df$Gene == sig_name, ]

pdf(file=file.path(dir.output, 'forestplot/PFS',
                    paste("INSPIRE_IOKIN", "pan_dicho.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logHR estimate",
                 label = "logHR")

dev.off()

################################################
## Forestplot: PFS Lung
################################################

df <- read.csv(file = file.path(dir.input, "sig_pfs.csv"))
df <- df[df$Gene == sig_name & df$Cancer_type == 'Lung', ]

pdf(file=file.path(dir.output, 'forestplot/PFS',
                    paste("INSPIRE_IOKIN", "lung_cont.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logHR estimate",
                 label = "logHR")

dev.off()


df <- read.csv(file = file.path(dir.input, "sig_pfs_dicho.csv"))
df <- df[df$Gene == sig_name & df$Cancer_type == 'Lung', ]

pdf(file=file.path(dir.output, 'forestplot/PFS',
                    paste("INSPIRE_IOKIN", "lung_dicho.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logHR estimate",
                 label = "logHR")

dev.off()


################################################
## Forestplot: Response
################################################
# pan-cancer
df <- read.csv(file = file.path(dir.input, "sig_logreg.csv"))
df <- df[df$Gene == sig_name, ]

pdf(file=file.path(dir.output, 'forestplot/response',
                    paste("INSPIRE_IOKIN", "pan.pdf", sep="_")),
     width = 10, height = 10)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logOR estimate",
                 label = "logOR")

dev.off()

# Melanoma
df <- read.csv(file = file.path(dir.input, "sig_logreg.csv"))
df <- df[df$Gene == sig_name &  df$Cancer_type == 'Melanoma', ]

pdf(file=file.path(dir.output, 'forestplot/response',
                    paste("INSPIRE_IOKIN", "melanoma.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logOR estimate",
                 label = "logOR")

dev.off()


# Lung
df <- read.csv(file = file.path(dir.input, "sig_logreg.csv"))
df <- df[df$Gene == sig_name &  df$Cancer_type == 'Lung', ]

pdf(file=file.path(dir.output, 'forestplot/response',
                    paste("INSPIRE_IOKIN", "lung.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logOR estimate",
                 label = "logOR")

dev.off()


# Kidney
df <- read.csv(file = file.path(dir.input, "sig_logreg.csv"))
df <- df[df$Gene == sig_name &  df$Cancer_type == 'Kidney', ]

pdf(file=file.path(dir.output, 'forestplot/response',
                    paste("INSPIRE_IOKIN", "kidney.pdf", sep="_")),
     width = 10, height = 8)

forestPlot(coef = df$Coef,
                 se = df$SE,
                 study  = substr(df$Study, 5, nchar(df$Study)),
                 pval = df$Pval,
                 n = df$N,
                 cancer.type = df$Cancer_type,
                 treatment = df$Treatment,
                 feature = unique(df$Gene),
                 xlab = "logOR estimate",
                 label = "logOR")

dev.off()

