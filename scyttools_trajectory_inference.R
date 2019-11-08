#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_trajectory_inference.R (-h | --help | --version)
scyttools_trajectory_inference.R RDS OUT

Description:   This script performs trajectory inference

Options:
--version       Show the current version.

Arguments:
RDS    Provide single cell experiment object Rdata file
OUT    Provide output file

" -> doc


args <- docopt(doc)

source("scyttools_functions.R")

library(monocle)

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

# for each supergroup, run monocle

cds_traj <- lapply(seq_along(unique(sce_glm_pca$supergroups)), function(supergroup){
  sce_supergroup <- sce_glm_pca[,sce_glm_pca$supergroups == supergroup]
  cds <- convertTo(sce_supergroup, type = "monocle",
                   row.fields = c(1:ncol(rowData(sce_supergroup))),
                   col.fields = c(1:ncol(colData(sce_supergroup))),
                   expressionFamily = negbinomial.size())
  
  featureData(cds)$gene_short_name <- rowData(sce_supergroup)$Symbol
  
  cds <- estimateSizeFactors(cds) 
  cds <- estimateDispersions(cds) 
  cds <- detectGenes(cds)
  
  cds <- setOrderingFilter(cds, row.names(fData(cds))[fData(cds)$dev <= 2000])
  cds <- reduceDimension(cds, method = "DDRTree")
  cds <- orderCells(cds)
  return(cds)
})

cds_col_data <- lapply(cds_traj, function(cds){
  col_data <- data.frame(cell_index = cds$cell_index,
                         monocle_pseudotime = cds$Pseudotime,
                         monocle_state = cds$State)
  return(col_data)
}) %>% bind_rows()

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  rownames_to_column("row_names") %>% 
  left_join(cds_col_data) %>%
  column_to_rownames("row_names") %>% 
  DataFrame()

DDRTree_matrix <- data.frame(barcode = rownames(reducedDim(sce_glm_pca, "UMAP"))) %>% 
  left_join(lapply(cds_traj, function(cds){
    t(monocle::reducedDimS(cds)) %>% 
      as.data.frame() %>% 
      rownames_to_column("barcode") %>% 
      return()}) %>%
    bind_rows()) %>% 
  column_to_rownames("barcode")

reducedDim(sce_glm_pca, "DDRTree") <- DDRTree_matrix

## write out new SingleCellExperiment Object
save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################