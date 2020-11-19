#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_integrate.R (-h | --help | --version)
scyttools_integrate.R OUT RDS RDS...

Description:   This script integrates multiple single cell objects

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

sce_list <- lapply(args$RDS,
                        function(sce_glm_pca){
                            load(sce_glm_pca)
                            sce_glm_pca <- sce
                            sce_glm_pca <- sce_glm_pca[order(rownames(rowData(sce_glm_pca))),]
                            rowData(sce_glm_pca) <- rowData(sce_glm_pca)[,c("ID", "Symbol")]
                            return(sce_glm_pca)
  }) %>%
  setNames(basename(dirname(args$RDS)))

rescaled <- multiBatchNorm(sce_list)

uncorrected <- do.call(cbind, rescaled)

sce <- devianceFeatureSelection(uncorrected,
  assay="counts",
  sorted=TRUE)

sce2<-sce[1:5000, ]
sce3<-GLMPCA(sce2,
  assay="counts",
  30)

glmpca_matrix <- as.matrix(metadata(sce3)$glmpca$factors)
batch <- factor(sce3$Sample %>%
    as.character() %>%
    dirname() %>%
    dirname() %>%
    basename())

mnn.out <- reducedMNN(glmpca_matrix,
  batch = batch,
  k=20)

sce_glm_pca <- sce
reducedDim(sce_glm_pca, "GLM_PCA") <- glmpca_matrix
reducedDim(sce_glm_pca, "corrected") <- mnn.out$corrected
sce_glm_pca$batch <- batch

sce_glm_pca <- runTSNE(sce_glm_pca, dimred="corrected")
sce_glm_pca <- runUMAP(sce_glm_pca, dimred="corrected")

g <- buildSNNGraph(sce_glm_pca, k=10, use.dimred = 'corrected')
clust <- igraph::cluster_walktrap(g)$membership

cell_data_set <- new_cell_data_set(expression_data = counts(sce_glm_pca),
                                   cell_metadata = colData(sce_glm_pca),
                                   gene_metadata = rowData(sce_glm_pca))

rownames(colData(cell_data_set)) <- paste(rownames(colData(cell_data_set)), colData(cell_data_set)$batch, sep = "_")

reducedDim(cell_data_set, "UMAP") <- reducedDim(sce_glm_pca, "UMAP")

cell_data_set <- cluster_cells(cell_data_set,
                               reduction_method = "UMAP",
                               partition_qval = 0.01)

colData(sce_glm_pca)$clust <- factor(clust)
colData(sce_glm_pca)$supergroups <- factor(cell_data_set@clusters[["UMAP"]]$partitions)
colData(sce_glm_pca)$monocle_clusters <- factor(cell_data_set@clusters[["UMAP"]]$clusters)

## write out new SingleCellExperiment Object

save(sce_glm_pca,
     file = args$OUT)