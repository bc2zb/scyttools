#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_integrate.R (-h | --help | --version)
scyttools_integrate.R RDS OUT

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

# load single cell data
load(args$RDS)

sce_list <- c("/data/capaldobj/CS027005_sc_rna_seq_cell_tag//SCAF1426_CT35-1-CTL_1w/trajectoried-sce.Rdata",
    "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1510_CT35-1-CTL_4w/trajectoried-sce.Rdata",
    "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1511_CT35-1-Ctrl/trajectoried-sce.Rdata")

sce_list <-lapply(sce_list,
                        function(sce_glm_pca){
                            load(sce_glm_pca)
                            reducedDim(sce_glm_pca, "GLM_PCA") <- as.matrix(reducedDim(sce_glm_pca, "GLM_PCA"))
                            sce_glm_pca <- sce_glm_pca[order(rownames(rowData(sce_glm_pca))),]
                            rowData(sce_glm_pca) <- rowData(sce_glm_pca)[,-4]
                            return(sce_glm_pca)
  }) %>% 
  setNames(basename(dirname(sce_list)))

dec_list <- lapply(sce_list, modelGeneVar)

combined.dec <- combineVar(dec_list[[1]], dec_list[[2]], dec_list[[3]])
chosen.hvgs <- combined.dec$bio > 0
summary(chosen.hvgs)

library(batchelor)

rescaled <- multiBatchNorm(sce_list) 

uncorrected <- cbind(rescaled[[1]], rescaled[[2]], rescaled[[3]])

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

scores <- c("/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1426_CT35-1-CTL_1w/geneset-scored-sce.Rdata",
  "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1510_CT35-1-CTL_4w/geneset-scored-sce.Rdata",
  "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1511_CT35-1-Ctrl/geneset-scored-sce.Rdata") %>% 
  lapply(function(sce_glm_pca){
    load(sce_glm_pca)
    scores <- colData(sce_glm_pca) %>% 
      as.data.frame() %>% 
      rownames_to_column("barcode") %>% 
      dplyr::select(barcode, NE_score, AR_score, cyclone_phase, monocle_pseudotime)
    return(scores)
  }) %>%
  setNames(c("SCAF1426_CT35-1-CTL_1w",
             "SCAF1510_CT35-1-CTL_4w",
             "SCAF1511_CT35-1-Ctrl")) %>% 
  bind_rows(.id = "batch")

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  rownames_to_column("rownames") %>% 
  dplyr::select(-c(clust, supergroups, monocle_clusters, cyclone_phase, monocle_pseudotime)) %>%
  left_join(scores, by = c("Barcode" = "barcode", "batch")) %>%
  column_to_rownames("rownames") %>% 
  DataFrame()


cell_data_set <- new_cell_data_set(expression_data = counts(sce_glm_pca),
                                   cell_metadata = colData(sce_glm_pca),
                                   gene_metadata = rowData(sce_glm_pca))

reducedDim(cell_data_set, "UMAP") <- reducedDim(sce_glm_pca, "UMAP")
cell_data_set <- cluster_cells(cell_data_set)
cell_data_set <- learn_graph(cell_data_set)
colData(sce_glm_pca)$supergroups <- factor(cell_data_set@clusters[["UMAP"]]$partitions)
colData(sce_glm_pca)$monocle_clusters <- factor(cell_data_set@clusters[["UMAP"]]$clusters)