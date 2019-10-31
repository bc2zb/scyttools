library(DropletUtils)
library(scater)
library(scran)
library(annotables)
library(tidyverse)
library(glmpca)
library(kohonen)
library(ggraph)
library(tidygraph)
library(igraph)
library(Seurat)
library(edgeR)
library(msigdbr)

set.seed(8675309)
source("scyttools_functions.R")

# zmad scoring cells

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


# load in signatures
neuro <- read_csv("~/lgcp/rnaseq/analysis-scripts/neuro_reference_vpca_loadings_grch37.csv")
ar_signature_weights <- read_csv("~/lgcp/rnaseq/analysis-scripts/ar_signature_weights_Mendiratta.csv")


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

# NE score ranking
logcounts(sce_glm_pca)[rowData(sce_glm_pca)$Symbol %in% neuro$Loading,] %>% 
  tidy_sparse_matrix() %>% 
  left_join(colData(sce_glm_pca) %>% 
              as.data.frame() %>% 
              select(Barcode, clust),
            by = c("column" = "Barcode")) %>% 
  left_join(neuro_joined, by = c("row" = "ID")) %>% 
  mutate(weighted_expression = value*PC1.v) %>% 
  group_by(clust, Loading) %>% 
  summarize(mean_weighted_expression = mean(weighted_expression, na.rm = T)) %>% 
  spread(clust,
         mean_weighted_expression) %>% 
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/NE-score-weighted-expression-values.csv")

# AR score ranking
logcounts(sce_glm_pca)[rownames(logcounts(sce_glm_pca)) %in% ar_signature_weights$`Ensemble ID`,] %>% 
  broom::tidy() %>% 
  left_join(colData(sce_glm_pca) %>% 
              as.data.frame() %>% 
              select(Barcode, clust),
            by = c("column" = "Barcode")) %>% 
  left_join(ar_signature_weights, by = c("row" = "Ensemble ID")) %>% 
  mutate(weighted_expression = value*`BinReg Coef`) %>% 
  group_by(clust, `Gene Symbol`) %>% 
  summarize(mean_weighted_expression = mean(weighted_expression, na.rm = T)) %>% 
  spread(clust,
         mean_weighted_expression) %>% 
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/AR-signature-score-weighted-expression-values.csv")
