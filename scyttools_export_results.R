#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_export_results.R (-h | --help | --version)
scyttools_export_results.R RDS OUT

Description:   This script performs trajectory inference

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

if(length(levels(sce_glm_pca$supergroups)) > 1){
  findMarkers(sce_glm_pca,
              sce_glm_pca$supergroups,
              direction = "up",
              lfc = 1) %>% 
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column, "ensgene") %>% 
    bind_rows(.id = "supergroup") %>% 
    left_join(grch38) %>% 
    dplyr::filter(Top <= 10) %>% 
  write_csv(paste0(args$OUT, "/partition-top-ranked-genes.csv"))
}

findMarkers(sce_glm_pca,
            sce_glm_pca$monocle_clusters,
            direction = "up",
            lfc = 1) %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "monocle cluster") %>% 
  left_join(grch38) %>% 
  dplyr::filter(Top <= 10) %>% 
write_csv(paste0(args$OUT, "/monocle-clusters-top-ranked-genes.csv"))

findMarkers(sce_glm_pca,
            sce_glm_pca$clust,
            direction = "up",
            lfc = 1) %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "monocle cluster") %>% 
  left_join(grch38) %>% 
  dplyr::filter(Top <= 10) %>% 
  write_csv(paste0(args$OUT, "/scyttools-clusters-top-ranked-genes.csv"))

reducedDim(sce_glm_pca, "UMAP") %>%
  as.data.frame() %>%
  rownames_to_column("Barcode") %>%
  setNames(c("Barcode", "UMAP-1", "UMAP-2")) %>%
  write_csv(paste0(args$OUT, "UMAP-projection.csv"))

colData(sce_glm_pca) %>%
  as.data.frame() %>%
  select(Barcode, monocle_clusters) %>% 
  write_csv(paste0(args$OUT, "monocle-clusters.csv"))

colData(sce_glm_pca) %>%
  as.data.frame() %>%
  select(Barcode, supergroups) %>% 
  write_csv(paste0(args$OUT, "partition-clusters.csv"))


colData(sce_glm_pca) %>%
  as.data.frame() %>%
  select(Barcode, clust) %>% 
  write_csv(paste0(args$OUT, "scyttools-clusters.csv"))

##########################################################################
############################     End code     ############################
##########################################################################
