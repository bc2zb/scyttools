#!/bin/bash

#Rscript scyttools.R --quality_control "/Volumes/Group09/CCB/Beshiri/Folders_old/Omics_data/CT35_10x_filtered_gbm/" "~/CT35-2-v1/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/CT35-2-v1/qc-sce.Rdata" "~/CT35-2-v1/dim-red-sce.Rdata"
#Rscript scyttools.R --cluster_cells "~/CT35-2-v1/dim-red-sce.Rdata" "~/CT35-2-v1/clustered-sce.Rdata"
#Rscript scyttools.R --trajectory_inference "~/CT35-2-v1/clustered-sce.Rdata" "~/CT35-2-v1/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/CT35-2-v1/trajectoried-sce.Rdata" "~/CT35-2-v1/geneset-scored-sce.Rdata"