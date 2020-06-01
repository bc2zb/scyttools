#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_cluster_cells.R (-h | --help | --version)
scyttools_cluster_cells.R RDS OUT

Description:   This script performs dimensionality reduction

Options:
--version       Show the current version.

Arguments:
RDS    Provide single cell experiment object Rdata file
OUT    Provide output file

" -> doc


args <- docopt(doc)

source("scyttools_functions.R")

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

# perform graph based clustering

g <- buildSNNGraph(sce_glm_pca, k=10, use.dimred = 'GLM_PCA')
clust <- igraph::cluster_walktrap(g)$membership

cell_data_set <- new_cell_data_set(expression_data = counts(sce_glm_pca),
                                   cell_metadata = colData(sce_glm_pca),
                                   gene_metadata = rowData(sce_glm_pca))

reducedDim(cell_data_set, "UMAP") <- reducedDim(sce_glm_pca, "UMAP")
cell_data_set <- cluster_cells(cell_data_set, partition_qval = 0.05)

colData(sce_glm_pca)$clust <- factor(clust)
colData(sce_glm_pca)$supergroups <- factor(cell_data_set@clusters[["UMAP"]]$partitions)
colData(sce_glm_pca)$monocle_clusters <- factor(cell_data_set@clusters[["UMAP"]]$clusters)

## write out new SingleCellExperiment Object

save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################