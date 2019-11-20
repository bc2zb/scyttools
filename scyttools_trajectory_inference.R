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

library(monocle)

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

# for each supergroup, run monocle

cds_traj <- lapply(seq_along(unique(sce_glm_pca$supergroups)), function(supergroup){
  sce_supergroup <- sce_glm_pca[,sce_glm_pca$supergroups == supergroup]
  
  hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
  cell_cycle_scores <- cyclone(sce_supergroup,
                               pairs = hs.pairs)
  
  cds <- convertTo(sce_supergroup, type = "monocle",
                   row.fields = c(1:ncol(rowData(sce_supergroup))),
                   col.fields = c(1:ncol(colData(sce_supergroup))),
                   expressionFamily = negbinomial.size())
  
  featureData(cds)$gene_short_name <- rowData(sce_supergroup)$Symbol
  
  cds <- estimateSizeFactors(cds) 
  cds <- estimateDispersions(cds) 
  cds <- detectGenes(cds)
  
  cds <- setOrderingFilter(cds, row.names(fData(cds))[fData(cds)$dev <= 2000])
  cds <- reduceDimension(cds, method = "DDRTree")
  cds <- orderCells(cds)
  
  cds$cyclone_phase <- cell_cycle_scores$phases
  cds$cyclone_G1_score <- cell_cycle_scores$scores$G1
  cds$cyclone_G2M_score <- cell_cycle_scores$scores$G2M
  cds$cyclone_S_score <- cell_cycle_scores$scores$S
  
  cds$cyclone_G1_score_normalized <- cell_cycle_scores$normalized.scores$G1
  cds$cyclone_G2M_score_normalized <- cell_cycle_scores$normalized.scores$G2M
  cds$cyclone_S_score_normalized <- cell_cycle_scores$normalized.scores$S
  
  cds_state_levels <- c(0, levels(cds$State))
  
  phase_state_counts <- pData(cds) %>%
    as.data.frame() %>%
    mutate(cyclone_phase = factor(cyclone_phase, levels = c("G1", "S", "G2M"))) %>% 
    group_by(cyclone_phase, State) %>%
    count() %>%
    ungroup() %>%
    mutate(State = factor(State, levels = cds_state_levels)) %>% 
    bind_rows(data.frame(cyclone_phase = c("G1", "S", "G2M"),                  # add in a fake state to ensure all phases are represented
                         State = factor(c(0,0,0), levels = cds_state_levels),
                         n = c(1,1,1))) %>% 
    spread(cyclone_phase, n, fill = 0) %>% 
    mutate(score = ((S+G2M)/(S+G2M+G1)),
           count = (S+G2M+G1),
           num_div = score*count) %>% 
    filter(State != 0) %>%                                                     # remove fake state
    mutate(State = factor(State, levels = levels(cds$State)))
  
  root_state <- phase_state_counts$State[which(phase_state_counts$num_div == max(phase_state_counts$num_div))]
  
  cds <- orderCells(cds, root_state = as.character(root_state))
  
  return(cds)
})

cds_col_data <- lapply(cds_traj, function(cds){
  col_data <- data.frame(cell_index = cds$cell_index,
                         monocle_pseudotime = cds$Pseudotime,
                         monocle_state = cds$State,
                         cyclone_phase = cds$cyclone_phase, 
                         cyclone_G1_score = cds$cyclone_G1_score,
                         cyclone_G2M_score = cds$cyclone_G2M_score,
                         cyclone_S_score = cds$cyclone_S_score,
                         cyclone_G1_score_normalized = cds$cyclone_G1_score_normalized,
                         cyclone_G2M_score_normalized = cds$cyclone_G2M_score_normalized,
                         cyclone_S_score_normalized = cds$cyclone_S_score_normalized)
  return(col_data)
}) %>% bind_rows()

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  rownames_to_column("row_names") %>% 
  left_join(cds_col_data) %>%
  column_to_rownames("row_names") %>% 
  DataFrame()

DDRTree_matrix <- data.frame(barcode = rownames(reducedDim(sce_glm_pca, "UMAP"))) %>% 
  left_join(lapply(cds_traj, function(cds){
    t(monocle::reducedDimS(cds)) %>% 
      as.data.frame() %>% 
      rownames_to_column("barcode") %>% 
      return()}) %>%
    bind_rows()) %>% 
  column_to_rownames("barcode")

reducedDim(sce_glm_pca, "DDRTree") <- DDRTree_matrix

## write out new SingleCellExperiment Object
save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################