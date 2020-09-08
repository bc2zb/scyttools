source("scyttools_functions.R")

sce_list <- c("/data/capaldobj/CS027005_sc_rna_seq_cell_tag//SCAF1427_CT35-2-CTL_1w/trajectoried-sce.Rdata",
    "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1512_CT35-2-CTL_4w/trajectoried-sce.Rdata",
    "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/SCAF1513_CT35-2-Ctrl/trajectoried-sce.Rdata")

load(sce_list[length(sce_list)])

seurat_list <-lapply(sce_list,
                        function(sce_glm_pca){
                            load(sce_glm_pca)
                            reducedDim(sce_glm_pca, "GLM_PCA") <- as.matrix(reducedDim(sce_glm_pca, "GLM_PCA"))
                            return(sce_glm_pca)
  }) %>% 
  setNames(basename(dirname(sce_list))) %>%
  lapply(as.Seurat) %>% 
  lapply(SCTransform)

features <- SelectIntegrationFeatures(object.list = seurat_list,
                                             nfeatures = 3000)

options(future.globals.maxSize = 4000 * 1024^2)

seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                         anchor.features = features, 
                                         verbose = TRUE)

seurat_list <- lapply(seurat_list, FindVariableFeatures, 
                             selection.method = "vst",
                             nfeatures = 2000,
                             verbose = TRUE)

anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                         normalization.method = "SCT",
                                         anchor.features = features,
                                         #reference = 1,
                                         dims = 1:30)

integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT",
                                   dims = 1:30)

integrated <- RunPCA(integrated,
    npcs = 30,
    verbose = TRUE,
    features = rowData(sce_glm_pca)$ID[1:5000])
integrated <- RunUMAP(integrated,
    reduction = "pca",
    dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.25)

DimPlot(integrated, reduction = "umap", label = TRUE) 
ggsave("umap-plot-ct35-2-colored-by-seurat-clusters.pdf") 
DimPlot(integrated, reduction = "umap", group.by = "Sample") + theme(legend.position = "none")
ggsave("umap-plot-ct35-2-colored-by-sample.pdf") 
save(integrated, file = "/data/capaldobj/CS027005_sc_rna_seq_cell_tag/seurat-integrated.Rdata") 