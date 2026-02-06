# -------------------------------------------------------------------------------
# Description:
# This script computes gene signature scores across multiple immune-oncology (IO) 
# transcriptomic datasets using a variety of established algorithms, including 
# singscore, GSVA, ssGSEA, weighted mean, and other signature-specific methods.
#
# For spatially informed (novel) signatures, samples are first scored using 
# rank-based single-sample gene set scoring and then projected into a 
# low-dimensional signature space defined by two predefined centroids derived 
# from spatial transcriptomics analyses.
#
# For each sample, Euclidean distances to both centroids are computed. The final 
# spatially informed signature score is defined as the distance to Centroid 2, 
# which represents a molecular phenotype associated with favorable 
# progression-free survival (PFS) and response to immune checkpoint blockade (ICB).
# Lower distances indicate greater similarity to the favorable spatial state.
#
# The script processes multiple preprocessed ICB datasets stored as 
# SummarizedExperiment objects and integrates signature scoring, centroid-based 
# distance computation (when applicable), and clinical outcome metadata.
#
# Required Inputs:
#   - SummarizedExperiment (.rda) files stored in `result/data`:
#       * Normalized gene expression matrices (e.g., log2(TPM + 1))
#       * Clinical outcome annotations (PFS time and event)
#   - Gene signature definitions in Ensembl or gene symbol format
#   - Signature metadata table (`signature_information.csv`) specifying scoring methods
#   - Two centroid profiles stored as .rds files (used for spatial signatures only)
#
# Parameters (parsed from file names):
#   - study_icb       – Unique study identifier (e.g., ICB_Gide)
#   - cancer_type     – Cancer type (e.g., Melanoma)
#   - treatment_type  – Immunotherapy treatment (e.g., PD-1/PD-L1)
#
# Analyses performed:
#   - Per-sample gene signature scoring using:
#       * singscore (simpleScore)
#       * GSVA
#       * ssGSEA
#       * Weighted Mean
#       * Signature-specific algorithms (e.g., COX-IS, IPS)
#   - Centroid-based distance scoring for spatially informed signatures
#
# Output:
#   - Signature score matrices saved as:
#       `result/score/<study>__<cancer>__<treatment>.rda`
#
# Notes:
#   - Centroid 2 corresponds to a spatial tumor–microenvironment state associated
#     with improved clinical outcomes under ICB.
#   - Signature gene identifiers must match those in the expression matrix.
#
# -------------------------------------------------------------------------------
#######################################################
# libraries
#######################################################

library(PredictioR)
library(dplyr)
library(survival)
library(data.table)
library(singscore)
library(MultiAssayExperiment)
library(GSVA)

source("workflow/external_validation_rnaseq/sigDistanceFunction.R")

#######################################################
# set up directory
#######################################################

dir.sig <- 'data/'
dir.input <- 'result/data'
dir.output <- 'result/score'

######################################################
# load centroids
######################################################
cent1 <- readRDS(file.path(dir.sig , 'centroid1_finalv2.rds'))
cent2 <- readRDS(file.path(dir.sig , 'Centroid2_finalv2.rds'))

####################################################################
## Load selected IO signature
####################################################################
files <- list.files(file.path(dir.input, 'symbol')) 

sig_file <- list.files(file.path('data', 'signature'))
sig_name <- substr(sig_file, 1, nchar(sig_file)-4)

signature_info <- read.csv(file.path('data', 'signature_information.csv')) 
int <- intersect(sig_name, signature_info$signature)

signature <- sapply(1:length(int), function(k){
  
  print(int[k])
  j <- which(sig_name == int[k])
  load(file.path('data', 'signature', sig_file[j]))
  sig

})

names(signature) <- int
###################################################################
## Compute signature score
###################################################################
## run signature score
for(k in 1:length(files)){ 

print(files[k])
load(file.path(dir.input, 'symbol', files[k]))  

# set up study name, cancer type, and treatment type
study_name <- substr(files[k], 5, nchar(files[k])-4)
study_icb <- strsplit(study_name, '__')[[1]][1]
cancer_type <- strsplit(study_name, '__')[[1]][2]
treatment_type <- strsplit(study_name, '__')[[1]][3]

# transform to log2(TPM+1)
annot <- rowData(dat_icb)
tpm <- assay(dat_icb)
tpm_updated <- (2 ^ tpm - 0.001)
expr <- log2(tpm_updated + 1)
if(sum(rownames(annot) == rownames(expr)) == nrow(expr)){

 rownames(expr) <- annot$gene_name

 
geneSig.score <- lapply(1:length(signature), function(i){ 
  
  print(paste(i , names(signature)[i], sep="/"))
  sig_name <- names(signature)[i]
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "GSVA"){
    
    geneSig <- geneSigGSVA(dat.icb = expr,
                           sig = signature[[i]],
                           sig.name = sig_name,
                           missing.perc = 0.5,
                           const.int = 1,
                           n.cutoff = 10,
                           sig.perc = 0.8,
                           study = study_icb,
                           gene.annot = "gene_name")
    
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Weighted Mean"){
    
    geneSig <- geneSigMean(dat.icb = expr,
                           sig = signature[[i]],
                           sig.name = sig_name,
                           missing.perc = 0.5,
                           const.int = 1,
                           n.cutoff = 10,
                           sig.perc = 0.8,
                           study = study_icb)
    
  }
  
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "ssGSEA"){
    
    geneSig <- geneSigssGSEA(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 1,
                             n.cutoff = 10,
                             sig.perc = 0.8,
                             study = study_icb)
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
    
  }
  
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita"){
    
    geneSig <- geneSigCOX_IS(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 1,
                             n.cutoff = 10,
                             sig.perc = 0.8,
                             study = study_icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong"){
    
    geneSig <- geneSigIPS(dat.icb = expr,
                          sig = signature[[i]],
                          sig.name = sig_name,
                          missing.perc = 0.5,
                          const.int = 1,
                          n.cutoff = 10,
                          study = study_icb)
    
  }
  
  #if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche"){
    
  #  geneSig <- geneSigPredictIO(dat.icb = expr,
  #                              sig = signature[[i]],
  #                              sig.name = sig_name,
  #                              missing.perc = 0.5,
  #                              const.int = 1,
  #                              n.cutoff = 10,
  #                              sig.perc = 0.8,
  #                              study = study_icb)
    
  #}
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo"){
    
    geneSig <- geneSigIPRES(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 1,
                            n.cutoff = 10,
                            sig.perc = 0.8,
                            study = study_icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du"){
    
    geneSig <- geneSigPassON(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 1,
                             n.cutoff = 10,
                             sig.perc = 0.8,
                             study = study_icb)
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen"){
    
    geneSig <- geneSigIPSOV(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 1,
                            n.cutoff = 10,
                            sig.perc = 0.8,
                            study = study_icb)
    
  }
  

  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "SpatialSignature_Marret"){
    
    geneSig <- geneSigNovel(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 1,
                            n.cutoff = 10,
                            centroid1 = cent1,
                            centroid2= cent2,
                            study = study_icb)
    
  }
  
  if(sum(!is.na(geneSig)) > 0){
    
    geneSig <- geneSig
    
  }     
  
  if(sum(!is.na(geneSig)) == 0){
    
    geneSig <- rep(NA, ncol(expr))
    
  }
  
  geneSig
  
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- names(signature)
remove <- which(is.na(rowSums(geneSig.score)))
if(length(remove) > 0){
  
  geneSig.score <- geneSig.score[-remove, ]
  
 }

save(geneSig.score, file=file.path(dir.output, paste(paste(study_icb, cancer_type, treatment_type, sep='__'), '.rda', sep="")))

}else{

 print('check annotation file')

}


}

