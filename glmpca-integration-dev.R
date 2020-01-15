source("scyttools_functions.R")

sce_file_list <- c("/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/SCAF1121_35_2A_CTL/qc-sce.Rdata",
                   "/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/SCAF1122_35_3A_CTL/qc-sce.Rdata",
                   "/Volumes/Group05/CCBB/Single-Cell-Bioinformatics-2019-October-03/SCAF1153_CT35_2_CTL/qc-sce.Rdata")

sce_list <- lapply(sce_file_list,
                   function(qc_sce_file){
                     load(qc_sce_file)
                     return(sce)
                   })

sce_list_rownames <- lapply(sce_list,
                            function(sce){
                              return(rownames(sce))
                            })
universe <- intersect(sce_list_rownames[[1]], intersect(sce_list_rownames[[2]], sce_list_rownames[[3]]))

sce_combined <- fastMNN(sce_list[[1]],
                        sce_list[[2]],
                        sce_list[[3]])
