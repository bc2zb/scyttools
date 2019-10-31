deg_all <- findMarkers(sce_glm_pca,
                   sce_glm_pca$clust,
                   direction = "up",
                   pval.type="all") %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "cluster_id") %>% 
  left_join(grch38)

deg_any <- findMarkers(sce_glm_pca,
                       sce_glm_pca$clust,
                       direction = "up",
                       pval.type="any") %>% 
  lapply(as.data.frame) %>% 
  lapply(rownames_to_column, "ensgene") %>% 
  bind_rows(.id = "cluster_id") %>% 
  left_join(grch38)

for(clust_id in unique(deg_all$cluster_id)){
  file <- paste0("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/", "ct-35-2-v1-10x-deg-all-cluster-id-", clust_id, ".csv")
  deg_all %>% 
    filter(cluster_id == clust_id) %>% 
    write_csv(file)
}

for(clust_id in unique(deg_any$cluster_id)){
  file <- paste0("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/", "ct-35-2-v1-10x-deg-any-cluster-id-", clust_id, ".csv")
  deg_any %>% 
    filter(cluster_id == clust_id) %>% 
    write_csv(file)
}

