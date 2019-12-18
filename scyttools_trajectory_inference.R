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

library(monocle3)

##########################################################################
############################ R code goes here ############################
##########################################################################

# load single cell data
load(args$RDS)

# add back cell cycle code
hs.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
cell_cycle_scores <- cyclone(sce_glm_pca,
                             pairs = hs.pairs)

colData(sce_glm_pca) <- colData(sce_glm_pca) %>% 
  as.data.frame() %>% 
  rownames_to_column("rownames") %>% 
  bind_cols(data.frame(cyclone_phase = cell_cycle_scores$phases)) %>% 
  column_to_rownames("rownames") %>% 
  DataFrame()

cell_data_set <- new_cell_data_set(expression_data = counts(sce_glm_pca),
                                   cell_metadata = colData(sce_glm_pca),
                                   gene_metadata = rowData(sce_glm_pca))

reducedDim(cell_data_set, "UMAP") <- reducedDim(sce_glm_pca, "UMAP")
cell_data_set <- cluster_cells(cell_data_set)
cell_data_set <- learn_graph(cell_data_set)

graph_info <- data.frame(vertex_name = igraph::V(principal_graph(cell_data_set)[["UMAP"]])$name,
                         vertex_degree = degree(principal_graph(cell_data_set)[["UMAP"]]),
                         component_membership = components(principal_graph(cell_data_set)[["UMAP"]])$membership)

# pull out neighborhood order number of vertices away from nodes
terminal_neighborhoods <- ego(cell_data_set@principal_graph$UMAP, order = 3, nodes = as.character((graph_info %>% filter(vertex_degree == 1))$vertex_name))
names(terminal_neighborhoods) <- as.character((graph_info %>% filter(vertex_degree == 1))$vertex_name)

if(length(terminal_neighborhoods) == 0){
  terminal_neighborhoods <- ego(cell_data_set@principal_graph$UMAP, order = 3, nodes = as.character((graph_info)$vertex_name))
  names(terminal_neighborhoods) <- as.character((graph_info)$vertex_name)
}

terminal_neighborhoods <- lapply(terminal_neighborhoods, function(neighborhood){return(data.frame(neighbors = as.vector(neighborhood)))}) %>%
  bind_rows(.id = "terminal_node")

root_nodes <- cell_data_set@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex %>% 
  as.data.frame() %>% 
  rownames_to_column("Barcode") %>% 
  left_join(colData(sce_glm_pca) %>% 
              as.data.frame()) %>% 
  group_by(cyclone_phase, V1, monocle_clusters) %>% 
  tally() %>% 
  ungroup() %>% 
  spread(cyclone_phase,
         n,
         fill = 0) %>% 
  mutate(vertex_name = paste0("Y_", V1)) %>% 
  left_join(graph_info) %>%
  left_join(terminal_neighborhoods %>% 
              mutate(vertex_name = paste0("Y_", neighbors))) %>% 
  filter(!is.na(terminal_node)) %>% 
  group_by(component_membership, terminal_node) %>% 
  summarize(num_dividing_cells = sum(S),
            num_G1 = sum(G1),
            num_G2M = sum(G2M),
            fraction_dividing = (num_dividing_cells + num_G2M)/(num_dividing_cells + num_G2M + num_G1)) %>% 
  ungroup() %>% 
  group_by(component_membership)

if(all(root_nodes$num_dividing_cells == 0)){
  root_nodes <- root_nodes %>% 
    summarise(root_node = terminal_node[which(fraction_dividing == max(fraction_dividing))])
}else{
  root_nodes <- root_nodes %>% 
    summarise(root_node = terminal_node[which(num_dividing_cells == max(num_dividing_cells))])
}
  
cell_data_set <- order_cells(cell_data_set, root_pr_nodes = root_nodes$root_node)
colData(sce_glm_pca)$monocle_pseudotime <- cell_data_set@principal_graph_aux$UMAP$pseudotime

## write out new SingleCellExperiment Object
save(sce_glm_pca,
          file = args$OUT)

##########################################################################
############################     End code     ############################
##########################################################################