source("scyttools_functions.R")

sce_list <- c("/data/capaldobj/CS027005_sc_rna_seq_cell_tag//SCAF1427_CT35-2-CTL_1w/trajectoried-sce.Rdata",
    "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1512_CT35-2-CTL_4w/trajectoried-sce.Rdata",
    "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1513_CT35-2-Ctrl/trajectoried-sce.Rdata")

sce_list <-lapply(sce_list,
                        function(sce_glm_pca){
                            load(sce_glm_pca)
                            reducedDim(sce_glm_pca, "GLM_PCA") <- as.matrix(reducedDim(sce_glm_pca, "GLM_PCA"))
                            sce_glm_pca <- sce_glm_pca[order(rownames(rowData(sce_glm_pca))),]
                            rowData(sce_glm_pca) <- rowData(sce_glm_pca)[,-4]
                            return(sce_glm_pca)
  }) %>% 
  setNames(basename(dirname(sce_list)))

dec_list <- lapply(sce_list, modelGeneVar)

combined.dec <- combineVar(dec_list[[1]], dec_list[[2]], dec_list[[3]])
chosen.hvgs <- combined.dec$bio > 0
summary(chosen.hvgs)

library(batchelor)

rescaled <- multiBatchNorm(sce_list) 

uncorrected <- cbind(rescaled[[1]], rescaled[[2]], rescaled[[3]])

sce <- devianceFeatureSelection(uncorrected,
  assay="counts",
  sorted=TRUE)
sce2<-sce[1:5000, ]
sce3<-GLMPCA(sce2,
  assay="counts",
  30)

glmpca_matrix <- as.matrix(metadata(sce3)$glmpca$factors)
batch <- factor(sce3$Sample %>%
    as.character() %>%
    dirname() %>% 
    dirname() %>% 
    basename())

mnn.out <- reducedMNN(glmpca_matrix,
  batch = batch,
  k=20)

sce_glm_pca <- sce
reducedDim(sce_glm_pca, "GLM_PCA") <- glmpca_matrix
reducedDim(sce_glm_pca, "corrected") <- mnn.out$corrected
sce_glm_pca$batch <- batch

snn.gr <- buildSNNGraph(sce_glm_pca, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn

sce_glm_pca <- runTSNE(sce_glm_pca, dimred="corrected")

plotTSNE(sce_glm_pca, colour_by="batch")

ggsave("mnn-corrected-tsne-plot.pdf")

reducedDim(sce_glm_pca, "TSNE") %>%
    as.data.frame() %>%
    bind_cols(data.frame(batch = sce_glm_pca$batch,
        cluster = factor(clusters.mnn))) %>%
    ggplot(aes(V1, V2, color = cluster)) +
    geom_point() +
    facet_grid(.~batch) +
    theme_bw()

ggsave("mnn-corrected-tsne-plot-facet-by-batch.pdf",
    width = 9,
    height = 4)

sce_glm_pca <- runUMAP(sce_glm_pca, dimred="corrected")

plotReducedDim(sce_glm_pca, "UMAP", colour_by="batch")

ggsave("mnn-corrected-umap-plot.pdf")

reducedDim(sce_glm_pca, "UMAP") %>%
    as.data.frame() %>%
    bind_cols(data.frame(batch = sce_glm_pca$batch,
        cluster = factor(clusters.mnn))) %>%
    ggplot(aes(V1, V2, color = cluster)) +
    geom_point() +
    facet_grid(.~batch) +
    theme_bw()

ggsave("mnn-corrected-umap-plot-facet-by-batch.pdf",
    width = 9,
    height = 4)

scores <- c("/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1512_CT35-2-CTL_4w/geneset-scored-sce.Rdata",
  "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1427_CT35-2-CTL_1w/geneset-scored-sce.Rdata",
  "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1513_CT35-2-Ctrl/geneset-scored-sce.Rdata") %>% 
  lapply(function(sce_glm_pca){
    load(sce_glm_pca)
    scores <- colData(sce_glm_pca) %>% 
      as.data.frame() %>% 
      rownames_to_column("barcode") %>% 
      dplyr::select(barcode, NE_score, AR_score, cyclone_phase, monocle_pseudotime)
    return(scores)
  }) %>%
  setNames(c("SCAF1512_CT35-2-CTL_4w",
             "SCAF1427_CT35-2-CTL_1w",
             "SCAF1513_CT35-2-Ctrl")) %>% 
  bind_rows(.id = "batch")

colData(sce_glm_pca) <- colData(sce_glm_pca) %>%
  as.data.frame() %>%
  rownames_to_column("rownames") %>% 
  dplyr::select(-c(clust, supergroups, monocle_clusters, cyclone_phase, monocle_pseudotime)) %>%
  left_join(scores, by = c("Barcode" = "barcode", "batch")) %>%
  column_to_rownames("rownames") %>% 
  DataFrame()


cell_data_set <- new_cell_data_set(expression_data = counts(sce_glm_pca),
                                   cell_metadata = colData(sce_glm_pca),
                                   gene_metadata = rowData(sce_glm_pca))

reducedDim(cell_data_set, "UMAP") <- reducedDim(sce_glm_pca, "UMAP")
cell_data_set <- cluster_cells(cell_data_set)
cell_data_set <- learn_graph(cell_data_set)
colData(sce_glm_pca)$supergroups <- factor(cell_data_set@clusters[["UMAP"]]$partitions)
colData(sce_glm_pca)$monocle_clusters <- factor(cell_data_set@clusters[["UMAP"]]$clusters)

save.image("ct35-2-integration.Rdata")