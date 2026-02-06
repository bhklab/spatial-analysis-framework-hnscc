# -------------------------------------------------------------------------------
# Description:
# This script evaluates the association between spatially informed gene signature
# metrics and clinical outcomes across multiple immunotherapy cohorts.
# The primary predictor is the Euclidean distance of each sample to a predefined
# reference centroid (centroid2), representing a spatial TME state associated with
# favorable progression-free survival (PFS) following immune checkpoint blockade.
#
# Signature distances are assumed to be pre-computed and stored in `.rda` files
# for each dataset.
#
# The script performs:
#   - Cox proportional hazards regression to assess association with PFS
#     (using distance to centroid2 as a continuous or dichotomized predictor)
#   - Logistic regression to assess association with objective response
#     (responder vs. non-responder)
#
# Required Inputs:
#   - Signature distance files (`.rda`) in `result/score/`, each containing:
#       * Distance to centroid2 per sample
#       * PFS time (`survival_time_pfs`) and event (`event_occurred_pfs`)
#       * Response status (`response`), when available
#
# Parameters:
#   - Dataset metadata (study, cancer type, treatment) parsed from file names
#   - Signature names derived from object row names
#
# Analyses performed:
#   - PFS association:
#       * Cox regression using centroid distance as a continuous predictor
#       * Cox regression using median-dichotomized distance (high vs. low)
#   - Response association:
#       * Logistic regression using centroid distance as a continuous predictor
#
# Output:
#   - `sig_pfs.rda` / `sig_pfs.csv` — Cox regression (continuous distance)
#   - `sig_pfs_dicho.rda` / `sig_pfs_dicho.csv` — Cox regression (dichotomized)
#   - `sig_logreg.rda` / `sig_logreg.csv` — Logistic regression (response)
#   - Outputs saved under `result/assoc/`
#
# Notes:
#   - Smaller distance to centroid2 indicates higher similarity to the
#     favorable spatial phenotype
#   - PFS is right-censored at 24 months
#   - Dichotomized models require ≥3 samples per group
#   - Response models require ≥1 sample per response category
# -------------------------------------------------------------------------------
############################################
## libraries
############################################
library(survival)
library(dplyr)
library(PredictioR)
library(data.table)
library(survcomp)
library(MultiAssayExperiment)

dir.input <- 'data/results/score'
dir.output <- 'data/results/'

############################################
## load signature score data
############################################
files <- list.files(dir.input)

geneSig.score <- lapply(1:length(files), function(k){
  load(file.path(dir.input, files[k]))
  geneSig.score
})

study_icb <- substr(files, 1, nchar(files) - 4)
names(geneSig.score) <- study_icb

#######################################
## Association with PFS (Cont)
#######################################
## load data
files <- list.files(file.path('data/results/data')) # 26 cohorts

## association continious
res_pfs <- lapply(1:length(files), function(k){ 

print(files[k])
load(file.path('data/results/data', files[k]))  

study_name <- substr(files[k], 1, nchar(files[k])-4)
study_icb <- strsplit(study_name, '__')[[1]][1]
cancer_type <- strsplit(study_name, '__')[[1]][2]
treatment_type <- strsplit(study_name, '__')[[1]][3]
var.pfs <- colData(dat_icb)$event_occurred_pfs

if(sum(!is.na(var.pfs)) > 0 & sum(var.pfs != 0, na.rm =TRUE) > 0 ){
 
 res <- lapply(1:nrow(geneSig.score[[k]]), function(i){

      geneSigSurvCont(dat.icb = dat_icb,
                      geneSig = geneSig.score[[k]][i, ],
                      time.censor = 24,
                      n.cutoff = 10,
                      study =  paste(study_icb, cancer_type, treatment_type, sep="__"),
                      surv.outcome = "PFS",
                      sig.name = rownames(geneSig.score[[k]])[i],
                      cancer.type = cancer_type,
                      treatment = treatment_type)

   })

 res <- do.call(rbind, res)

}else{


    res <-  data.frame(Outcome = 'PFS',
                       Gene = NA,
                       Study = paste(study_icb, cancer_type, treatment_type, sep="__"),
                       Coef = NA,
                       SE = NA,
                       N = NA,
                       Pval = NA,
                       Cancer_type = NA,
                       Treatment = NA 
                      #FDR = NA
                       )

}
   
#if(nrow(res) != 1 & sum(!is.na(res$Coef)) > 0){
#res <- res[!is.na(res$Coef), ]
#res$FDR <- p.adjust(res$Pval, method = 'BH')
#}

res

})

res_pfs <- do.call(rbind, res_pfs)
res_pfs <- res_pfs[!is.na(res_pfs$Coef), ]

save(res_pfs, file = file.path(dir.output, 'assoc', 'sig_pfs.rda'))
write.csv(res_pfs, file = file.path(dir.output, 'assoc', 'sig_pfs.csv'), row.names=FALSE)

#######################################
## Association with PFS (Dicho)
#######################################
## load data
files <- list.files(file.path('data/results/data'))

## association dicho (High vs Low)
res_pfs <- lapply(1:length(files), function(k){ 

print(files[k])
load(file.path('data/results/data', files[k]))  

study_name <- substr(files[k], 1, nchar(files[k])-4)
study_icb <- strsplit(study_name, '__')[[1]][1]
cancer_type <- strsplit(study_name, '__')[[1]][2]
treatment_type <- strsplit(study_name, '__')[[1]][3]
clin <- colData(dat_icb)


if(sum(!is.na(clin$event_occurred_pfs)) > 0 &  sum(clin$event_occurred_pfs != 0, na.rm =TRUE) > 0 ){

   res <- lapply(1:nrow(geneSig.score[[k]]), function(i){

  
   cox <- survDicho(status = clin$event_occurred_pfs,
                   time = clin$survival_time_pfs, 
                   time.censor = 24,
                   var = geneSig.score[[k]][i, ],
                   n0.cutoff = 3,
                   n1.cutoff = 3,
                   method = "median", 
                   var.type = FALSE)

         data.frame(Outcome = "PFS", 
                    Gene = rownames(geneSig.score[[k]])[i],
                    Study = paste(study_icb, cancer_type, treatment_type, sep="__"), 
                    Coef = round(cox["HR"], 3), 
                    SE = round(cox["SE"],3), 
                    N = cox["N"], 
                    Pval = cox["Pval"], 
                    Cancer_type = cancer_type,        
                    Treatment = treatment_type)     

})

res <- do.call(rbind, res)

}else{


    res <-  data.frame(Outcome = 'PFS',
                       Gene = NA,
                       Study = paste(study_icb, cancer_type, treatment_type, sep="__"),
                       Coef = NA,
                       SE = NA,
                       N = NA,
                       Pval = NA,
                       Cancer_type = NA,
                       Treatment = NA
                      #FDR = NA
                       )

}
   
#if(nrow(res) != 1 & sum(!is.na(res$Coef)) > 0){
#res <- res[!is.na(res$Coef), ]
#res$FDR <- p.adjust(res$Pval, method = 'BH')
#}

res

})

res_pfs <- do.call(rbind, res_pfs)
res_pfs <- res_pfs[!is.na(res_pfs$Coef), ]

save(res_pfs, file = file.path(dir.output, 'assoc', 'sig_pfs_dicho.rda'))
write.csv(res_pfs, file = file.path(dir.output, 'assoc', 'sig_pfs_dicho.csv'), row.names=FALSE)

#######################################
## Association with response (R/NR)
#######################################
## load data
files <- list.files(file.path('data/results/data'))

## association response
res_logreg <- lapply(1:length(files), function(k){ 

print(files[k])
load(file.path('data/results/data', files[k]))  

study_name <- substr(files[k], 1, nchar(files[k])-4)
study_icb <- strsplit(study_name, '__')[[1]][1]
cancer_type <- strsplit(study_name, '__')[[1]][2]
treatment_type <- strsplit(study_name, '__')[[1]][3]
updated_sig <- geneSig.score[[k]]
#updated_sig <- updated_sig[rowSums(updated_sig == 0) != ncol(updated_sig), ]

res <- lapply(1:nrow(updated_sig), function(i){

   geneSigLogReg(dat.icb = dat_icb,
                 geneSig = updated_sig[i, ],
                 n.cutoff = 10,
                 study =  paste(study_icb, cancer_type, treatment_type, sep="__"),
                 sig.name = rownames(updated_sig)[i],
                 n0.cutoff = 1, 
                 n1.cutoff =  1,
                 cancer.type = cancer_type,
                 treatment = treatment_type)

})

res <- do.call(rbind, res)
#res <- res[!is.na(res$Coef), ]
#res$FDR <- p.adjust(res$Pval, method = 'BH')
res

})

res_logreg <- do.call(rbind, res_logreg)

save(res_logreg, file = file.path(dir.output, 'assoc', 'sig_logreg.rda'))
write.csv(res_logreg, file = file.path(dir.output, 'assoc', 'sig_logreg.csv'), row.names=FALSE)
