#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_dimensionality_reduction.R (-h | --help | --version)
scyttools_dimensionality_reduction.R RDS OUT

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

sce_d <- sce[1:5000,]
sce_d <- GLMPCA(sce_d, assay="counts", 30)

reducedDim(sce, "GLM_PCA") <- reducedDim(sce_d, "GLMPCA")

sce_glm_pca <- runUMAP(sce, dimred = "GLM_PCA", pca = 30)

## write out new SingleCellExperiment Object

save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################