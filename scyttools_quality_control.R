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

location_tidy <- rowData(sce) %>% 
  as.data.frame() %>% 
  left_join(grch38 %>% 
              select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
  distinct()
is.mito <- which(location_tidy$chr == "MT")

if(length(is.mito) == 0){
  location_tidy <- rowData(sce) %>% 
    as.data.frame() %>% 
    left_join(grcm38 %>% 
                select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
    distinct()
  is.mito <- which(location_tidy$chr == "MT")
}


sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))

# remove cells with log-library sizes that are more than 3 MADs below the median
qc.lib <- isOutlier(log(sce$total), nmads=3, type="lower")
# remove cells where the log-transformed number of expressed genes is 3 MADs below the median
qc.nexprs <- isOutlier(log(sce$detected), nmads=3,
                       type="lower")
# remove cells where the number of mito genes is 3 MADs above the median
qc.mito <- isOutlier(sce$subsets_Mito_percent, nmads=3, type="higher")

discard <- qc.lib | qc.nexprs | qc.mito

DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs), MitoProp=sum(qc.mito), Total=sum(discard)) 

# Retain only high-quality cells in the SingleCellExperiment.
sce <- sce[,!discard]

## normalization

cl<-scran::quickCluster(sce)
sce<-scran::computeSumFactors(sce,clusters=cl)
sce <- scater::logNormCounts(sce)

## feature selection

colnames(sce) <- colData(sce)$Barcode

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