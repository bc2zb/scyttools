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

set.seed(8675309)

# source functions
source("scyttools_functions.R")

# load single cell data, this should be passed in from command line
sce <- read10xCounts("/Volumes/Group09/CCB/Beshiri/Folders_old/CT35/'Omics_data/single_cell_RNAseq/lineage-tracing-May-2018/CT35_10x_filtered_gbm")

# should be loading reference from commmand line
location_tidy <- rowData(sce) %>% 
  as.data.frame() %>% 
  left_join(grch38 %>% 
              select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
  distinct()
is.mito <- which(location_tidy$chr == "MT")

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

sce_d <- sce[ranked_genes$dev <= 2000,]

filtered_counts <- counts(sce_d)
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0,]

glmpca_poi_30 <- glmpca(as.matrix(filtered_counts),
                        30,
                        fam = "poi")

reducedDim(sce, "GLM_PCA") <- as.matrix(glmpca_poi_30$factors)

sce_glm_pca <- runUMAP(sce, use_dimred = "GLM_PCA", pca = 30)

rowData(sce_glm_pca) <- cbind(rowData(sce_glm_pca),
                              ranked_genes)

rowData(sce_glm_pca) <- rowData(sce_glm_pca) %>% 
  as.data.frame() %>% 
  rownames_to_column("row_names") %>% 
  left_join(cbind(glmpca_poi_30$loadings, glmpca_poi_30$coefX) %>% 
              rownames_to_column("row_names") %>% 
              rename("coefX" = "V1")) %>%
  column_to_rownames("row_names") %>% 
  DataFrame()

# perform graph based clustering

g <- buildSNNGraph(sce_glm_pca, k=10, use.dimred = 'GLM_PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

colData(sce_glm_pca)$clust <- factor(clust)
colData(sce_glm_pca)$NE_score <- barcode_neuro_score
colData(sce_glm_pca)$AR_score <- barcode_ar_signature_score

# perform som-kmeans-signature-scoring

full_som_grid <- somgrid(xdim = 30, ydim=30, topo="rectangular")

full_som_model <- som(reducedDim(sce, "GLM_PCA"), 
                      grid=full_som_grid, 
                      rlen=10000, 
                      alpha=c(0.05,0.01), 
                      keep.data = TRUE )

som_grid <- somgrid(xdim = 2, ydim=2, topo="rectangular")
# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available

som_results <- lapply(seq_along(1:(ncol(reducedDim(sce, "GLM_PCA"))-1)),
                      function(col_index){
                        dim_x <- reducedDim(sce, "GLM_PCA")[,col_index]
                        cols_to_test_start <- col_index + 1
                        cols_to_test_end <- ncol(reducedDim(sce, "GLM_PCA"))
                        colnames_index <- cols_to_test_start:cols_to_test_end
                        dims_y <- reducedDim(sce, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                        dims_y <- as.data.frame(dims_y)
                        colnames(dims_y) <- colnames(reducedDim(sce, "GLM_PCA"))[colnames_index]
                        
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

kmeans_results <- lapply(seq_along(1:(ncol(reducedDim(sce, "GLM_PCA"))-1)),
                         function(col_index){
                           dim_x <- reducedDim(sce, "GLM_PCA")[,col_index]
                           cols_to_test_start <- col_index + 1
                           cols_to_test_end <- ncol(reducedDim(sce, "GLM_PCA"))
                           colnames_index <- cols_to_test_start:cols_to_test_end
                           dims_y <- reducedDim(sce, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                           dims_y <- as.data.frame(dims_y)
                           colnames(dims_y) <- colnames(reducedDim(sce, "GLM_PCA"))[colnames_index]
                           
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

split_results <- lapply(seq_along(1:(ncol(reducedDim(sce, "GLM_PCA"))-1)),
                        function(col_index){
                          dim_x <- reducedDim(sce, "GLM_PCA")[,col_index]
                          cols_to_test_start <- col_index + 1
                          cols_to_test_end <- ncol(reducedDim(sce, "GLM_PCA"))
                          colnames_index <- cols_to_test_start:cols_to_test_end
                          dims_y <- reducedDim(sce, "GLM_PCA")[,cols_to_test_start:cols_to_test_end]
                          dims_y <- as.data.frame(dims_y)
                          colnames(dims_y) <- colnames(reducedDim(sce, "GLM_PCA"))[colnames_index]
                          
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

save(sce_glm_pca, file = "/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/ct-35-2-v1-sce-object.Rdata")

reducedDim(sce_glm_pca, "UMAP") %>%
  as.data.frame() %>%
  rownames_to_column("Barcode") %>%
  setNames(c("Barcode", "UMAP-1", "UMAP-2")) %>%
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/ct35-2-v1-UMAP-projection.csv")

colData(sce_glm_pca) %>%
  as.data.frame() %>%
  select(Barcode, clust) %>% 
  write_csv("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/ct35-2-v1-glm-pca-graph-clusters.csv")