# -------------------------------------------------------------------------------
# Description:
# This script performs correlation and meta-correlation analysis of gene 
# signature scores (or distances) across multiple immunotherapy datasets. 
#
# The analysis assumes signature scores or distance metrics are precomputed 
# and stored as `.rda` files (one per dataset) in the `result/score/` folder.
#
# Main analyses include:
#   - Pearson correlation of signature scores within each dataset
#   - Meta-analysis of pairwise correlations across studies (e.g., pan-cancer, Lung or Melanoma)
#   - Integration of results into a meta-correlation matrix
#   - Hierarchical clustering to identify co-regulated or co-occurring signatures
#   - Heatmap generation with cluster annotations and signature types
#
# Required Inputs:
#   - Signature score `.rda` files per study in `result/score/`
#   - Signature metadata (`signature_information.csv`)
#
# Parameters:
#   - Study names are parsed from file names
#   - Analysis subset can be filtered by cancer type (e.g., Lung or Melanoma)
#
# Outputs:
#   - Per-study Pearson correlation matrices
#   - Meta-correlations across selected studies
#   - Meta-correlation matrix
#   - Signature clustering results
#   - Annotated heatmap PDF
#
# Notes:
#   - Meta-correlation requires ≥3 datasets for each pair of signatures
#   - Hierarchical clustering is used to define signature clusters
#   - Heatmap includes annotations for signature type and assigned cluster
# -------------------------------------------------------------------------------
############################################
## libraries
############################################
library(dplyr)
library(Hmisc)
library(meta)
library(NbClust)
library(colorRamp2)
library(PredictioR)
library(data.table)
library(ComplexHeatmap)

dir.input <- 'data/results/score'
dir.output <- 'data/results/cor'

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

# signature information
sig.info <- read.csv(file.path('data/results/sig', 'signature_information.csv')) 

############################################################
## Pearson Correlation analysis
############################################################

for(k in 1:length(geneSig.score)){

print(study_icb[k])
fit <- rcorr(t(geneSig.score[[k]]))
cor.val <- list ("r" = fit$r,
                 "n" = fit$n,
                 "p" = fit$P)

names(cor.val) <- rep(study_icb[k], 3)
save(cor.val, file = file.path(dir.output, 'cor', paste(study_icb[k], "sig_pcor.rda", sep = "_")))

}

################################################
## Merge correlation data
################################################
cor.files <- list.files(file.path(dir.output, 'cor'))

cor <- lapply(1:length(cor.files), function(i){

  print(i)
  load(file.path(dir.output, 'cor', cor.files[i]))
  study <- names(cor.val)[1]
  r <- cor.val[[1]]
  n <- cor.val[[2]]
  p <- cor.val[[3]]

  cor.res <- lapply(1:nrow(r), function(k){

    print(k)
    data.frame(study = study,
               geneSig1 = rownames(r)[k],
               geneSig2 = colnames(r),
               r = r[k, ],
               n = n[k, ],
               p = p[k, ])

  })

  cor.res <- do.call(rbind, cor.res)
  rownames(cor.res) <- NULL
  cor.res

})

cor <- do.call(rbind, cor)
################################################
## Correlation integration
################################################
group <- 'Kidney'
pcor <- cor
# Select only those that contain "Melanoma" & "Lung"
# pcor <- cor[grepl("Melanoma", cor$study), ]
# pcor <- cor[grepl("Lung", cor$study), ]
 pcor <- cor[grepl("Kidney", cor$study), ]

geneSig1 <- unique(pcor$geneSig1)

meta.cor <- lapply(1:length(geneSig1), function(i){

  print(i)
  df <- pcor[pcor$geneSig1 == geneSig1[i], ]
  geneSig2 <- unique(df$geneSig2)

 res <- lapply(1:length(geneSig2), function(j){

      print(j)
      sub.df <- df[df$geneSig2 == geneSig2[j], ]

      if(geneSig2[j] == geneSig1[i]){

        meta.res <- data.frame(geneSig1 = geneSig1[i],
                               geneSig2 = geneSig2[j],
                               r = 1,
                               se = NA,
                               pval = NA,
                               I2 = NA)
       }

     if(sum(!is.na(sub.df$p)) > 0 & nrow(sub.df) >= 3){

      fit <- metacor(sub.df$r, sub.df$n, sm = "cor",
                     control=list(stepadj=0.5, maxiter=1000))
      meta.res <- data.frame(geneSig1 = geneSig1[i],
                             geneSig2 = geneSig2[j],
                             r = fit$TE.random,
                             se = fit$seTE.random,
                             pval = fit$pval.random,
                             I2 = fit$I2)
    }

    if(nrow(sub.df) < 3){

      meta.res <- data.frame(geneSig1 = geneSig1[i],
                             geneSig2 = geneSig2[j],
                             r = NA,
                             se = NA,
                             pval = NA,
                             I2 = NA)


    }

    meta.res

  })

 do.call(rbind, res)
})

meta.cor <- do.call(rbind, meta.cor)
save(meta.cor, file=file.path(dir.output, group, "metaCor.RData"))

#########################################################
## Create a matrix of meta correlation
#########################################################
load(file.path(dir.output, group, "metaCor.RData"))

meta.cor <- meta.cor[!is.na(meta.cor$r), ]
geneSig1 <- unique(meta.cor$geneSig1)
geneSig2 <- unique(meta.cor$geneSig1)

metaCor <- lapply(1:length(geneSig1), function(i){

  print(i)
  df <- meta.cor[meta.cor$geneSig1 == geneSig1[i], ]

  cor.res <- sapply(1:length(geneSig2), function(j){

    print(j)

    if(sum(df$geneSig2 == geneSig2[j]) > 0){
       r <- unique(df[df$geneSig2 == geneSig2[j], "r"])
       if(r > 1){ r <- 1 }
       if(r < (-1)){ r <- (-1)}
     }else{ r <- NA }

    r

  })

  cor.res

})

metaCor <- do.call(cbind, metaCor)
rownames(metaCor) <- geneSig1
colnames(metaCor) <- geneSig2

save(metaCor, file=file.path(dir.output, group, "metaCor_matrix.RData"))

###################################################
## Hierarchical clustering
###################################################
load(file.path(dir.output, group, "metaCor_matrix.RData"))

metaCor <- metaCor[, colnames(metaCor) %in% sig.info$signature]
metaCor <- metaCor[rownames(metaCor) %in% sig.info$signature, ]

## find optimal number of clusters
res <- NbClust(metaCor, distance = "euclidean", min.nc=1, max.nc=15, 
             method = "complete", index = "ch")

res$Best.nc[1]

# Compute the distance matrix
dist_matrix <- dist(metaCor)

# Perform hierarchical clustering
row.clus<-hclust(dist(metaCor, method = "euclidean"), method = "complete")
col.clus<-hclust(dist(t(metaCor), method = "euclidean"), method = "complete")

# Plot the dendrogram (optional)

pdf(file=file.path(dir.output, group, "ht_sig.pdf"),
     width = 10, height = 8)
plot(row.clus, labels = rownames(metaCor), 
     main = "Dendrogram", xlab = "", sub = "")
dev.off()

# Extract clusters by specifying the number of clusters
clusters_k <- cutree(row.clus, k = res$Best.nc[1])  # Change 'k' to your desired number of clusters 

# Extract clusters by specifying the height
clusters_h <- cutree(col.clus, k = res$Best.nc[1])  # Change 'h' to your desired height

# View clusters
cluster <- lapply(1:length(clusters_h), function(k){
  data.frame( signature_name = names(clusters_h[k]),
              clust = clusters_h[k])
})

cluster <- do.call(rbind, cluster)
write.csv(cluster, file = file.path(dir.output, group, "cluster_cor.csv"), row.names = F )

############################################
## Heatmap for meta correlation data
############################################

cluster.cor <- read.csv(file.path(dir.output, group, "cluster_cor.csv"))
load(file.path(dir.output, group, "metaCor_matrix.RData"))

int <- intersect(sig.info$signature, cluster.cor$signature_name)

df <- metaCor
df <- df[order(rownames(df)), ]
df <- df[, order(colnames(df))]
df <- df[rownames(df) %in% int, ]
df <- df[, colnames(df) %in% int]

sig.info <- sig.info[sig.info$signature %in% rownames(df), ]
sig.info <- sig.info[order(sig.info$signature), ]
sig.info <- sig.info[sig.info$signature %in% int, ]

cluster.cor <- cluster.cor[cluster.cor$signature %in% sig.info$signature, ]
cluster.cor <- cluster.cor[order(cluster.cor$signature_name), ]
cluster.cor <- cluster.cor[cluster.cor$signature %in% int, ]

rownames(df) <- sig.info$signature
colnames(df) <- sig.info$signature
cluster.cor$signature_name <- sig.info$signature

col = list("Signature type" = c("resistance" = "#855C75FF", "sensitive" = "#526A83FF", "novel" = "#853820ff"),
           "Cluster" = c( "1" = "#78847FFF" , "2" = "#B1AF53FF", "3" = "#8491BEFF"
                       #   "4" = "#64894DFF", "5" = "#9C6755FF","6" = "#8491BEFF",
                       #   "7" = "#D2C396FF",  "8" = "#735c18ff"
                       ) 
           ) 

#Create the heatmap annotation
ha <- HeatmapAnnotation(
  "Signature type" = sig.info$association,
  "Cluster" = cluster.cor$clust,
   col = col
)

col_hclust <- hclust(dist(t(df), method = "euclidean"), method = "complete")
row_hclust <- hclust(dist(t(df), method = "euclidean"), method = "complete")

symmetric_matrix[upper.tri(symmetric_matrix)] <- NA

# Combine the heatmap and the annotation

pdf(file=file.path(dir.output, group, "ht_cor_sig.pdf"),
     width = 6, height = 5)

ht_data <- Heatmap(df, name = "Cor",
                   cluster_rows = as.dendrogram(row_hclust),
                   cluster_columns = as.dendrogram(col_hclust),
                   top_annotation = ha,
                   #right_annotation = row_annot,
                   show_row_names = FALSE, 
                   show_column_names = FALSE,
                   colorRamp2(c(-1, 0, 1), c("#4685A0FF", "white", "#864568FF")))

ComplexHeatmap::draw(ht_data,
     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right")

dev.off()


