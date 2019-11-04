#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools_cluster_cells.R (-h | --help | --version)
scyttools_cluster_cells.R DIR

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

# perform graph based clustering

g <- buildSNNGraph(sce_glm_pca, k=10, use.dimred = 'GLM_PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

colData(sce_glm_pca)$clust <- factor(clust)

# perform som-kmeans-signature-scoring

full_som_grid <- somgrid(xdim = 30, ydim=30, topo="rectangular")

full_som_model <- som(reducedDim(sce_glm_pca, "GLM_PCA"), 
                      grid=full_som_grid, 
                      rlen=10000, 
                      alpha=c(0.05,0.01), 
                      keep.data = TRUE )

som_grid <- somgrid(xdim = 2, ydim=2, topo="rectangular")
# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available

som_results <- lapply(seq_along(1:(ncol(reducedDim(sce_glm_pca, "GLM_PCA"))-1)),
                      function(col_index){
                        dim_x <- reducedDim(sce_glm_pca, "GLM_PCA")[,col_index]
                        cols_to_test_start <- col_index + 1
                        cols_to_test_end <- ncol(reducedDim(sce_glm_pca, "GLM_PCA"))
                        colnames_index <- cols_to_test_start:cols_to_test_end
                        dims_y <- reducedDim(sce_glm_pca, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                        dims_y <- as.data.frame(dims_y)
                        colnames(dims_y) <- colnames(reducedDim(sce_glm_pca, "GLM_PCA"))[colnames_index]
                        
                        subs_som_results <- lapply(dims_y,
                                                   function(dim_y){
                                                     som_data <- data.frame(dim_x = dim_x,
                                                                            dim_y = dim_y) %>% 
                                                       as.matrix()
                                                     som_model <- som(som_data, 
                                                                      grid=som_grid, 
                                                                      rlen=10000, 
                                                                      alpha=c(0.05,0.01), 
                                                                      keep.data = TRUE )
                                                     return(som_model$unit.classif)
                                                   })
                        return(subs_som_results)
                      }) %>% 
  lapply(bind_cols) %>%
  lapply(rownames_to_column, "cell_index")
names(som_results) <- paste0("som_dim", 1:length(som_results))

som_results <- som_results %>% 
  bind_rows(.id = "dim_x") %>% 
  gather(dim_y,
         som_code,
         -c(dim_x, cell_index)) %>% 
  filter(!is.na(som_code)) %>% 
  mutate(som_code = factor(som_code)) %>% 
  unite(dim_x_y,
        dim_x,
        dim_y) %>% 
  spread(dim_x_y,
         som_code)

kmeans_results <- lapply(seq_along(1:(ncol(reducedDim(sce_glm_pca, "GLM_PCA"))-1)),
                         function(col_index){
                           dim_x <- reducedDim(sce_glm_pca, "GLM_PCA")[,col_index]
                           cols_to_test_start <- col_index + 1
                           cols_to_test_end <- ncol(reducedDim(sce_glm_pca, "GLM_PCA"))
                           colnames_index <- cols_to_test_start:cols_to_test_end
                           dims_y <- reducedDim(sce_glm_pca, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                           dims_y <- as.data.frame(dims_y)
                           colnames(dims_y) <- colnames(reducedDim(sce_glm_pca, "GLM_PCA"))[colnames_index]
                           
                           subs_kmeans_results <- lapply(dims_y,
                                                         function(dim_y){
                                                           kmeans_data <- data.frame(dim_x = dim_x,
                                                                                     dim_y = dim_y) %>% 
                                                             as.matrix()
                                                           kmeans_clust <- kmeans(kmeans_data, centers = 4, iter.max = 1000)
                                                           return(kmeans_clust$cluster)
                                                         })
                           return(subs_kmeans_results)
                         }) %>% 
  lapply(bind_cols) %>%
  lapply(rownames_to_column, "cell_index")
names(kmeans_results) <- paste0("kmeans_dim", 1:length(kmeans_results))

kmeans_results <- kmeans_results %>% 
  bind_rows(.id = "dim_x") %>% 
  gather(dim_y,
         som_code,
         -c(dim_x, cell_index)) %>% 
  filter(!is.na(som_code)) %>% 
  mutate(som_code = factor(som_code)) %>% 
  unite(dim_x_y,
        dim_x,
        dim_y) %>% 
  spread(dim_x_y,
         som_code)

split_results <- lapply(seq_along(1:(ncol(reducedDim(sce_glm_pca, "GLM_PCA"))-1)),
                        function(col_index){
                          dim_x <- reducedDim(sce_glm_pca, "GLM_PCA")[,col_index]
                          cols_to_test_start <- col_index + 1
                          cols_to_test_end <- ncol(reducedDim(sce_glm_pca, "GLM_PCA"))
                          colnames_index <- cols_to_test_start:cols_to_test_end
                          dims_y <- reducedDim(sce_glm_pca, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                          dims_y <- as.data.frame(dims_y)
                          colnames(dims_y) <- colnames(reducedDim(sce_glm_pca, "GLM_PCA"))[colnames_index]
                          
                          subs_kmeans_results <- lapply(dims_y,
                                                        function(dim_y){
                                                          kmeans_data <- data.frame(dim_x = dim_x,
                                                                                    dim_y = dim_y) %>% 
                                                            rownames_to_column("barcode") %>% 
                                                            gather(dim,
                                                                   position,
                                                                   starts_with("dim")) %>% 
                                                            group_by(barcode) %>% 
                                                            summarize(distance = (sum(position^2))^(1/2))
                                                          kmeans_clust <- kmeans(kmeans_data$distance, centers = 2, iter.max = 1000)
                                                          if(kmeans_clust$centers[1,1] > kmeans_clust$centers[2,1]){
                                                            clusters <- data.frame(clust_id = kmeans_clust$cluster) %>% 
                                                              mutate(clust_id = if_else(clust_id == 1, 1, 0))
                                                          }else{
                                                            clusters <- data.frame(clust_id = kmeans_clust$cluster) %>% 
                                                              mutate(clust_id = if_else(clust_id == 2, 1, 0))
                                                          }
                                                          return(clusters$clust_id)
                                                        })
                          return(subs_kmeans_results)
                        }) %>% 
  lapply(bind_cols) %>%
  lapply(rownames_to_column, "cell_index")
names(split_results) <- paste0("split_dim", 1:length(split_results))
split_results <- split_results %>% 
  bind_rows(.id = "dim_x") %>% 
  gather(dim_y,
         som_code,
         -c(dim_x, cell_index)) %>% 
  filter(!is.na(som_code)) %>% 
  mutate(som_code = factor(som_code)) %>% 
  unite(dim_x_y,
        dim_x,
        dim_y) %>% 
  spread(dim_x_y,
         som_code)

colData(sce_glm_pca)$cell_index <- as.character(c(1:nrow(colData(sce_glm_pca))))
colData(sce_glm_pca)$full_som_node <- full_som_model$unit.classif

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  rownames_to_column("row_names") %>% 
  left_join(som_results) %>%
  left_join(kmeans_results) %>%
  left_join(split_results) %>%
  column_to_rownames("row_names") %>% 
  DataFrame()

## write out new SingleCellExperiment Object

save(sce_glm_pca,
     paste(RESULTS_DIR, "scyttools_sce_cluster_cells.Rdata"))

##########################################################################
############################     End code     ############################
##########################################################################