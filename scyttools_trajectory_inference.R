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

library(monocle3)

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

cell_data_set <- new_cell_data_set(expression_data = counts(sce_glm_pca),
                                   cell_metadata = colData(sce_glm_pca),
                                   gene_metadata = rowData(sce_glm_pca))

reducedDim(cell_data_set, "UMAP") <- reducedDim(sce_glm_pca, "UMAP")
cell_data_set <- cluster_cells(cell_data_set)
cell_data_set <- learn_graph(cell_data_set)

igraph::V(principal_graph(cell_data_set)[["UMAP"]])$name

cell_data_set@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex



## write out new SingleCellExperiment Object
save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################