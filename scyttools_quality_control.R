#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_quality_control.R (-h | --help | --version)
scyttools_quality_control.R DIR OUT

Description:   This script performs quality control on 10x genomics scRNA gene barcode matrix data from the cellranger pipeline

Options:
--version       Show the current version.

Arguments:
DIR    Provide directory where scyttools.args.Rdata file is located
OUT    Provide output file

" -> doc


args <- docopt(doc)

source("scyttools_functions.R")

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
sce <- read10xCounts(args$DIR)
sce <- scDblFinder(sce,
                   verbose = F)

sce <- sce[,sce$scDblFinder.class != "doublet"]
# should be loading reference from commmand line
# identify genes on mitochondria to filter out cells high in MT genes

qc_metrics <- perCellQCMetrics(sce)

## normalization

cl<-scran::quickCluster(sce)
sce<-scran::computeSumFactors(sce,clusters=cl)
sce <- scater::logNormCounts(sce)

## feature selection

colnames(sce) <- colData(sce)$Barcode
colData(sce) <- cbind(colData(sce), qc_metrics)
ranked_genes <- rank_all_genes(sce, "total")

## append ranked genes to singleCellExperiment object

rowData(sce) <- cbind(rowData(sce),
                              ranked_genes)

## write out new SingleCellExperiment Object

save(sce,
     file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################