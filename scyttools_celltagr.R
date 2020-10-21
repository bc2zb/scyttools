#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_celltagr.R (-h | --help | --version)
scyttools_celltagr.R OUT WHITELIST INPUT

Description:   This script performs celltagr analysis

Arguments:

WHITELIST   celltag whitelist
INPUT       input csv file containing rows of single cell object and bam file
OUT         Provide output directory

" -> doc

args <- docopt(doc)

library(CellTagR)

input_table <- read_csv(args$INPUT)

dir.create(args$OUT, showWarnings = F, recursive = T)
setwd(args$OUT)
dir.create("combined_bams/", showWarnings = F, recursive = T)
for(i in 1:nrow(input_table)){
    sub_wd <- paste("sample", i, sep = "_")
    dir.create(sub_wd, showWarnings = F, recursive = T)
    setwd(sub_wd)
    load(input_table$sce[i])
    write_tsv(data.frame(barcodes = sce$Barcode),
        path = "barcodes.tsv",
        col_names = F)
    command <- paste0("samtools view -f 12 -b ", shQuote(input_table$bam[i]), " > unmapped_read_pairs.bam")
    system(command)
    command <- paste0("samtools view -b  ", shQuote(input_table$bam[i]), " EGFP > egfp_mapped_read_pairs.bam")
    system(command)
    command <- paste0("samtools merge celltags.bam unmapped_read_pairs.bam egfp_mapped_read_pairs.bam")
    system(command)
    command <- paste0("cp celltags.bam ../combined_bams/", sub_wd, ".bam")
    system(command)
    rm(sce)
    setwd(args$OUT)
}

celltag <- CellTagObject("celltagged",
                         fastq.bam.directory = "combined_bams")

celltag <- CellTagExtraction(celltag.obj = celltag,
                             celltag.version = "v1")

celltag_whitelist <- args$WHITELIST

barcode_files <- paste0("sample_", 1:nrow(input_table), "/barcodes.tsv")

Barcode.Aggregate(barcode_files, "./barcodes_all.tsv")

celltag <- CellTagMatrixCount(celltag, barcodes.file = "./barcodes_all.tsv")

dir.create("star_collapse/", showWarnings = F, recursive = T)

celltag <- CellTagDataForCollapsing(celltag, "./collapsing.txt")

for(i in 1:nrow(input_table)){
    command <- paste0("~/starcode/starcode -s --print-clusters _Sample-", i, ".txt > star_collapse/collapsing_result_Sample-", i, ".txt")
    system(command)
}

celltag <- CellTagDataPostCollapsing(celltag,
                                     collapsed.rslt.file = list.files("star_collapse/",
                                                                      full.names = T))

celltag <- SingleCellDataBinatization(celltag, 2)

celltag <- SingleCellDataWhitelist(celltag, celltag_whitelist)

celltag <- MetricBasedFiltering(celltag, 20, comparison = "less")
celltag <- MetricBasedFiltering(celltag, 2, comparison = "greater")

celltag <- JaccardAnalysis(celltag)
celltag <- CloneCalling(celltag, 0.7)

save(celltag, file = "celltag-object.Rdata")
