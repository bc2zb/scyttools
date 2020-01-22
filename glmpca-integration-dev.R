source("scyttools_functions.R")
library(harmony)
sce_file_list <- c("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/SCAF1120_35_1A_CTL/qc-sce.Rdata",
                   "/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/SCAF1152_CT35_1_CTL/qc-sce.Rdata")

sce_list <- lapply(sce_file_list,
                   function(qc_sce_file){
                     load(qc_sce_file)
                     #sce <- SingleCellExperiment(assays = list(counts = counts(sce)))
                     split_file <- str_split(qc_sce_file, "/")
                     run_id_index <- str_which(split_file[[1]], "qc-sce.Rdata")
                     run_id <- split_file[[1]][run_id_index - 1]
                     colnames(sce) <- paste0(run_id,
                                             ".", colnames(sce))
                     colData(sce) <- colData(sce)[,1:2]
                     rowData(sce) <- rowData(sce)[,1:3]
                     return(sce)
                   })

sce_combined <- do.call(cbind, sce_list)

location_tidy <- data.frame(ID = rownames(sce_combined)) %>% 
  left_join(grch38 %>% 
              select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
  distinct()
is.mito <- which(location_tidy$chr == "MT")

if(length(is.mito) == 0){
  location_tidy <- rowData(sce_combined) %>% 
    as.data.frame() %>% 
    left_join(grcm38 %>% 
                select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
    distinct()
  is.mito <- which(location_tidy$chr == "MT")
}

# should be passing in multicore param from command line
# QC steps, could probably be a separate script, then QC'd data can be passed around to each of the subroutines
sce_combined <- calculateQCMetrics(sce_combined, feature_controls=list(Mito=is.mito), BPPARAM = MulticoreParam(7))
colData(sce_combined) <- cbind(colData(sce_combined), data.frame(run_id = str_remove(rownames(colData(sce_combined)), "\\.[A-Z]*-1$")))

ranked_genes <- rank_all_genes(sce_combined, "total_counts")

rowData(sce_combined) <- cbind(rowData(sce_combined),
                               ranked_genes)

sce_d <- sce_combined[rowData(sce_combined)$dev <= 2000,]

filtered_counts <- counts(sce_d)
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0,]

glmpca_poi_30 <- glmpca(as.matrix(filtered_counts),
                        30,
                        fam = "poi",
                        penalty = 1)

reducedDim(sce_combined, "GLM_PCA") <- as.matrix(glmpca_poi_30$factors)
reducedDim(sce_combined, "GLM_PCA_Corrected") <- as.matrix(HarmonyMatrix(as.matrix(glmpca_poi_30$factors),
                                                                         colData(sce_combined)$run_id,
                                                                         do_pca = F,
                                                                         max.iter.harmony = 100))

sce_glm_pca <- runUMAP(sce_combined, use_dimred = "GLM_PCA_Corrected", pca = 30)

