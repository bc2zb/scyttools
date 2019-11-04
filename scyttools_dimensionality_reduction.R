#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_dimensionality_reduction.R (-h | --help | --version)
scyttools_dimensionality_reduction.R DIR

Description:   This script performs dimensionality reduction

Options:
--version       Show the current version.

Arguments:
DIR    Provide directory where scyttools.args.Rdata file is located

" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

load(paste(ARGS_DIR, "scyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("scyttools_functions.R")

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

sce_d <- sce[rowData(sce)$dev <= 2000,]

filtered_counts <- counts(sce_d)
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0,]

glmpca_poi_30 <- glmpca(as.matrix(filtered_counts),
                        30,
                        fam = "poi")

reducedDim(sce, "GLM_PCA") <- as.matrix(glmpca_poi_30$factors)

sce_glm_pca <- runUMAP(sce, use_dimred = "GLM_PCA", pca = 30)

rowData(sce_glm_pca) <- rowData(sce_glm_pca) %>% 
  as.data.frame() %>% 
  rownames_to_column("row_names") %>% 
  left_join(cbind(glmpca_poi_30$loadings, glmpca_poi_30$coefX) %>% 
              rownames_to_column("row_names") %>% 
              rename("coefX" = "V1")) %>%
  column_to_rownames("row_names") %>% 
  DataFrame()

## write out new SingleCellExperiment Object

save(sce_glm_pca,
     paste(RESULTS_DIR, "scyttools_sce_dimensionality_reduction.Rdata"))

##########################################################################
############################     End code     ############################
##########################################################################