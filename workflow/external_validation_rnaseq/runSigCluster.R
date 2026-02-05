################################################################################
# Description:
#   Cluster immuno-oncology (IO) gene-expression signatures based on shared-gene
#   overlap, visualize signature similarity in PCA space, and annotate clusters
#   using KEGG pathway enrichment (Enrichr).
#
# What this script does
#   1) Loads all signature .rda files from: data/signature/
#        - each file must load an object named `sig`
#        - `sig` is either:
#            a) a data.frame with columns: signature_name, gene_name, OR
#            b) a vector/list of gene symbols
#   2) Filters signatures to those listed in: data/signature_information.csv
#   3) Builds a long-format table of all (signature_name, gene_name)
#   4) Computes a pairwise overlap matrix (signature x signature):
#        overlap[i, j] = number of shared genes between signature i and j
#      and saves it for reuse
#   5) Performs PCA on the overlap matrix (scaled), then runs affinity propagation
#      clustering on the first two PCs (PC1/PC2) to obtain signature clusters
#   6) Produces two outputs for visualization:
#        - PCA scatter plot with signature labels
#        - affinity propagation heatmap of clusters/similarity
#   7) Exports a signature-to-cluster mapping table as CSV
#   8) For each cluster with >1 signature:
#        - takes the union of genes across signatures in the cluster
#        - runs KEGG enrichment via Enrichr (KEGG_2016)
#        - saves the top 25 pathways per cluster as CSV
#
# Inputs
#   data/signature/                      : folder of signature .rda files
#     - each file loads an object `sig` (see format above)
#   data/signature_information.csv       : metadata table; must contain column:
#     - signature (used to intersect with signature filenames)
#
# Key parameters / assumptions
#   - Excludes CIN signatures by default:
#       "CIN25_Carter.rda", "CIN70_Carter.rda"
#   - Clustering is performed on PC1/PC2 coordinates only (not all PCs)
#   - Enrichment database is set to: KEGG_2016 (Enrichr)
#   - Gene identifiers are assumed to be human gene symbols (for Enrichr)
#
# Outputs (written under result/cluster/)
#   overlap.rda                         : signature x signature overlap matrix
#   IO_sig_gene.pdf                     : PCA scatter plot + labels
#   IO_sig_gene_heatmap.pdf             : affinity propagation clustering heatmap
#   clusters_IO_genes.csv               : two columns: signature, cluster
#   IO_cluster<k>.csv                   : KEGG enrichment results (top 25) per
#                                        cluster k (clusters with >1 signature)
################################################################################
###########################################
# Load libraries
###########################################
library(stringr)
library(apcluster)
library(RColorBrewer)
library(ggfortify)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(enrichR)

###########################################
## Signature information
###########################################

dir_in <- 'data'
dir_out <- 'result/cluster'

###########################################################
## Load IO data and signatures
###########################################################
sig_file <- list.files(file.path(dir_in, 'signature'))
sig_name <- substr(sig_file, 1, nchar(sig_file)-4)

signature_info <- read.csv(file.path(dir_in, 'signature_information.csv')) # 53 signatures
int <- intersect(sig_name, signature_info$signature)
sig_file <- sig_file[which(sig_name  %in% int)]
sig_name <- sig_name[which(sig_name  %in% int)]

signature <- sapply(1:length(sig_file), function(k){

  load(file.path(dir_in, 'signature', sig_file[k]))
  sig
})
names(signature) <- sig_name

#######################################################################################
## Identify common genes
#######################################################################################
dir_GeneSig <- file.path(dir_in, 'signature')
GeneSig_list <- list.files(dir_GeneSig)
GeneSig_list <- GeneSig_list[order(GeneSig_list)]
GeneSig_list <- GeneSig_list[!GeneSig_list %in% c("CIN25_Carter.rda", "CIN70_Carter.rda") ]

data <- lapply(1:length(GeneSig_list), function(k){

  load(file.path(dir_GeneSig, GeneSig_list[k]))
  
  if(class(sig) == 'data.frame'){
   df <- data.frame(signature_name = sig$signature_name,
             gene_name = sig$gene_name)
  }else{

   df <- data.frame(signature_name = substr(GeneSig_list[k], 1, nchar(GeneSig_list[k])-4),
             gene_name = unlist(sig))

  }

  df
 
})

data <- do.call(rbind, data)

# number of genes in common across signatures
signature <- sort( unique( data$signature_name ) )

overlap <- matrix( nrow = length( signature ) , ncol = length( signature ) , 0 )
colnames( overlap ) <- rownames( overlap ) <- signature

for( i in 1:length( signature ) ){

  print(i)

  for( j in 1:length( signature ) ){

    s1 <- data[ data$signature_name %in% signature[i] , ]$gene_name
    s2 <- data[ data$signature_name %in% signature[j] , ]$gene_name

    int <- intersect( s1 , s2 )

    overlap[ i , j ] <- length( int )

  }
}

save(overlap, file=file.path(dir_out, 'overlap.rda'))

##########################################################
## IO signatures
##########################################################
# PCA analysis
load(file=file.path(dir_out, 'overlap.rda'))

df <- overlap
pca <- prcomp( df , scale = TRUE)
var.res <- pca$sdev^2 / sum(pca$sdev^2)
x1 <- pca$x[ , 1:2]

apres <- apcluster(negDistMat(r=2), x1)

## Visualization plot
pdf(file = file.path(dir_out, "IO_sig_gene_final.pdf"), 
     width = 8, height = 5)

par(
  cex.axis = 1.1,   # size of axis tick labels
  cex.lab  = 1.3   # size of axis titles
)

plot( apres , x1 , xlab = "PC1 (21%)" , ylab = "PC2 (10%)" )
#text( x1, row.names(x1) , cex = 0.3 , pos = 4 , col = "black" )

dev.off()


## Visualization plot
#pdf(file = file.path(dir_out, "IO_sig_gene_heatmap.pdf"), 
#     width = 10, height = 10)

#heatmap(apres)
#text( x1, row.names(x1) , cex = 0.4 , pos = 4 , col = "black" )

#dev.off()

## save the results
cl <- apres@clusters
cluster <- NULL

for( i in 1:length( cl ) ){

  c <- unlist( cl[ i ] )
  cluster <- rbind( cluster , cbind( names( c ) , i ) )

}

cluster <- as.data.frame( cluster )
colnames(cluster ) <- c( "signature" , "cluster" )

cluster$signature <- as.character( cluster$signature )
cluster$cluster <- as.character( cluster$cluster )

write.csv( cluster , file = file.path(dir_out, "clusters_IO_genes.csv") )

############################################
## KEGG pathway analysis
############################################
setEnrichrSite("Enrichr") # Human genes
dbs <- c( "KEGG_2016" )

## IO signatures
cluster <- read.csv(file.path(dir_out, "clusters_IO_genes.csv"))
GeneSig_list <- GeneSig_list[order(GeneSig_list)]
GeneSig_list <- GeneSig_list[!GeneSig_list %in% c("CIN25_Carter.rda", "CIN70_Carter.rda") ]
signature_name <- substr(GeneSig_list, 1, nchar(GeneSig_list) -4)

clust <- sort( names( table( cluster$cluster )[ table( cluster$cluster ) > 1 ] ) )

int <- cluster$signature
sig <- lapply(1:length(int), function(k){
  
  print(k)
  load(file.path(dir_GeneSig, GeneSig_list[which(signature_name == int[k])]))
  if(class(sig) == 'data.frame'){

    df <- sig[, c("signature_name", "gene_name")]

  }else{

    df <- data.frame(signature_name = substr(GeneSig_list[k], 1, nchar(GeneSig_list[k])-4),
             gene_name = unlist(sig))

  }
  
 df

})

data <- do.call(rbind, sig)

for( i in 1:length( clust ) ){

  print(i)
  sig <- cluster[ cluster$cluster %in% clust[ i ] , ]$signature

  gene <- sort( unique( data[ data$signature_name %in% sig , ]$gene_name ) )
  enriched <- enrichr( gene, dbs)

  kegg <- cbind( clust[ i ] , enriched[[1]][ 1:25 , 1:7 ] )
  colnames(kegg) <- c("cluster", "pathway", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value",
                      "Old.Adjusted.P.value", "Odds.Ratio")
  write.csv(kegg, file= file.path(dir_out, paste(paste("IO_cluster", clust[i], sep=""), ".csv", sep="")),
             row.names = FALSE)

  }


