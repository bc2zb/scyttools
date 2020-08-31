#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_geneset_scoring.R (-h | --help | --version)
scyttools_geneset_scoring.R RDS OUT

Description:   This script performs trajectory inference

Options:
--version       Show the current version.

Arguments:
RDS    Provide single cell experiment object Rdata file
OUT    Provide output file

" -> doc


args <- docopt(doc)

source("scyttools_functions.R")

library(data.table)

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

# calculate mean absolute deviation for each gene across all samples
median_absolute_deviations <- logcounts(sce_glm_pca) %>% 
  tidy_sparse_matrix() %>%
  group_by(row) %>%
  summarise(median_abs_dev = mad(value),
            median = median(value))

# calculate modified Z-score by multipling the difference of the cpm_tmm and mean for each gene times 0.6745 (just because), and dividing by the mean absolute deviation
modified_z_score <- logcounts(sce_glm_pca) %>% 
  tidy_sparse_matrix() %>%
  left_join(median_absolute_deviations) %>%
  mutate(mod_z = (0.6745*(value - median))/median_abs_dev) %>% 
  select(row, column, mod_z)

fake_cell <- as.character(ncol(sce_glm_pca) + 1)

zmad_sparse <- data.frame(row = names(logcounts(sce_glm_pca)[,1])) %>% 
  left_join(modified_z_score) %>% 
  mutate(row = factor(row, levels = names(logcounts(sce_glm_pca)[,1])),
         column = if_else(is.na(column), fake_cell, column),
         mod_z = if_else(is.na(mod_z), 0, mod_z),
         column = factor(column, levels = c(colnames(logcounts(sce_glm_pca)), fake_cell))) %>% 
  arrange(column) %>% 
  tidytext::cast_sparse(row = row,
                        column = column,
                        value = mod_z)

assay(sce_glm_pca, "zmad") <- zmad_sparse[,colnames(zmad_sparse) != fake_cell]

geneset_table <- msigdbr::msigdbr() %>% 
  select(gs_name, entrez_gene) %>% 
  left_join(grch38 %>% 
              distinct(ensgene, entrez),
            by = c("entrez_gene" = "entrez")) %>% 
  setDT()

dt_zmad <- zmad_sparse %>% 
  tidy_sparse_matrix() %>% 
  setDT()

zmad_scores <- lapply(unique(geneset_table$gs_name), function(geneset){
  barcode_scores <- dt_zmad[row %in% geneset_table$ensgene[geneset_table$gs_name == geneset], mean(value, na.rm = T), by = column]  
  return(barcode_scores)
})

names(zmad_scores) <- unique(geneset_table$gs_name)
zmad_scores <- zmad_scores %>% 
  bind_rows(.id = "gs_name") %>% 
  spread(gs_name,
         V1)

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  rownames_to_column("row_names") %>% 
  left_join(zmad_scores,
            by = c("Barcode" = "column")) %>%
  column_to_rownames("row_names") %>% 
  DataFrame()

# load in signatures
neuro <- read_csv("/data/capaldobj/nf-core/lgcp/rnaseq/analysis-scripts/neuro_reference_vpca_loadings_grch37.csv")
ar_signature_weights <- read_csv("/data/capaldobj/nf-core/lgcp/rnaseq/analysis-scripts/ar_signature_weights_Mendiratta.csv")


barcode_ar_signature_score <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i],
                       ensgene = names(logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,i])) %>% 
    left_join(ar_signature_weights, by = c("ensgene" = "Ensemble ID")) %>% 
    transmute(score = logcounts*`BinReg Coef`)
  
  barcode_ar_signature_score[i] <- sum(scored$score)
}

neuro_joined <- neuro %>% 
  left_join(as.data.frame(rowData(sce_glm_pca)), by = c("Loading" = "Symbol"))

barcode_neuro_score <- vector(length = ncol(sce_glm_pca))
for(i in 1:ncol(sce_glm_pca)){
  scored <- data.frame(logcounts = logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,i],
                       ensgene = names(logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,i])) %>% 
    left_join(neuro_joined, by = c("ensgene" = "ID")) %>% 
    transmute(score = logcounts*PC1.v)
  
  barcode_neuro_score[i] <- sum(scored$score)
}

sce_glm_pca$NE_score <- barcode_neuro_score
sce_glm_pca$AR_score <- barcode_ar_signature_score

## write out new SingleCellExperiment Object
save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################