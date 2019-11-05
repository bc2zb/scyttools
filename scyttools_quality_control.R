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

# should be loading reference from commmand line
# identify genes on mitochondria to filter out cells high in MT genes
location_tidy <- rowData(sce) %>% 
  as.data.frame() %>% 
  left_join(grch38 %>% 
              select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
  distinct()
is.mito <- which(location_tidy$chr == "MT")

# should be passing in multicore param from command line
# QC steps, could probably be a separate script, then QC'd data can be passed around to each of the subroutines
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=is.mito), BPPARAM = MulticoreParam(7))

# remove cells with log-library sizes that are more than 3 MADs below the median
qc.lib <- isOutlier(log(sce$total_counts), nmads=3, type="lower")
# remove cells where the log-transformed number of expressed genes is 3 MADs below the median
qc.nexprs <- isOutlier(log(sce$total_features_by_counts), nmads=3,
                       type="lower")
# remove cells where the number of mito genes is 3 MADs above the median
qc.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")

discard <- qc.lib | qc.nexprs | qc.mito

# Retain only high-quality cells in the SingleCellExperiment.
sce <- sce[,!discard]

## normalization

cl<-scran::quickCluster(sce)
sce<-scran::computeSumFactors(sce,clusters=cl)
sce <- scater::normalize(sce)

## feature selection

colnames(sce) <- colData(sce)$Barcode
ranked_genes <- rank_all_genes(sce, "total_counts")

## append ranked genes to singleCellExperiment object

rowData(sce) <- cbind(rowData(sce),
                              ranked_genes)

## write out new SingleCellExperiment Object

save(sce,
     file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################