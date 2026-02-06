# ------------------------------------------------------------
# Spatially Informed Gene Signature Scoring via Centroid Distance
#
# This script defines helper functions to compute a spatially
# informed gene expression signature from bulk or spatial RNA-seq
# data. Signature activity is first quantified using rank-based
# single-sample scoring, after which samples are projected into a
# low-dimensional signature space defined by reference centroids.
#
# For each sample, Euclidean distances to two predefined centroids
# (e.g., representing distinct spatial tumor microenvironment
# states) are computed. The final signature score is defined as the
# Euclidean distance to a selected reference centroid (centroid 2),
# which in the discovery cohort corresponds to a spatial state
# associated with favorable progression-free survival following
# immune checkpoint blockade (ICB).
#
# Lower distances to centroid 2 indicate greater similarity to the
# favorable spatial TME state and are interpreted as higher
# predicted sensitivity to ICB.
#
# Functions:
#   - compute_distances(): Computes Euclidean distances from a
#     sample to two centroids in signature space.
#   - geneSigNovel(): Computes rank-based signature scores,
#     projects samples onto centroid space, and returns the
#     distance to the reference centroid (centroid 2) as the
#     spatially informed signature score per sample.
#
# Input data formats:
#   - SummarizedExperiment
#   - MultiAssayExperiment
#   - matrix or data.frame (genes × samples)
#
# Output:
#   - Numeric vector of centroid-distance–based signature scores
#     (one value per sample)
#
# Dependencies:
#   - singscore (rankGenes, simpleScore)
#   - SummarizedExperiment / MultiAssayExperiment
# ------------------------------------------------------------
############################################
## define function
############################################

compute_distances<-function(id, centroid1, centroid2){
  dist1<-sqrt(sum((id-centroid1)^2))
  dist2<-sqrt(sum((id-centroid2)^2))
  return(c(dist1,dist2))
}

geneSigNovel <- function (dat.icb, sig, sig.name, missing.perc = 0.50, const.int = 1, 
                          n.cutoff, study, centroid1, centroid2){
  
  if (!class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", 
                             "data.frame", "matrix")) {
    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))
  }
  
  if (class(dat.icb) == "MultiAssayExperiment") {
    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
  }
  
  if (class(dat.icb) == "SummarizedExperiment") {
    dat_expr <- assay(dat.icb)
  }
  
  if (class(dat.icb) %in% c("data.frame", "matrix")) {
    dat_expr <- dat.icb
  }
  
  data <- dat_expr
  #remove <- rem(data, missing.perc, const.int)
  
   # if (length(remove)) {
   #  data <- data[-remove, ]
   # }
  
  # scale.data <- scalefun(data)
  geneSig <- NULL
  
  if ( ncol(data) >= n.cutoff ) {
    
    
      # Step 1 ---- compute signature score
      ranked_expr <- rankGenes(data)
      sigScore_ssGSEA <- lapply(1:length(sig), function(j){
      
      sigScore <-simpleScore(ranked_expr, upSet = sig[[j]]) 
      # extract signature score
      sigScore$TotalScore
      
    })
    
    sigScore_ssGSEA <- do.call(cbind, sigScore_ssGSEA)
    
    if(ncol(sigScore_ssGSEA) != length(sig)){
      
      print("not enough signatures were computed")
      
    }else{
      
      rownames(sigScore_ssGSEA) <- colnames(data)
      colnames(sigScore_ssGSEA) <- names(sig)
      
      # Order data 
      centroid1 <- centroid1[, order(colnames(centroid1))]
      centroid2 <- centroid2[, order(colnames(centroid2))]
      sigScore_ssGSEA <- sigScore_ssGSEA[, order(colnames(sigScore_ssGSEA))]
      
      # Step 2 ---- compute distance
      
      dist <- t(apply(sigScore_ssGSEA, 1, compute_distances, centroid1, centroid2))

      # Focus on distance to centroid2, which in the training set is associated with more favorable PFS outcomes following immuno
      geneSig <- dist[,2]

    }
    
  }else{
       
    print("not enough samples")
    
  }
  
  return(geneSig)
  
}  
  
