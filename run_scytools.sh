#!/bin/bash

#Rscript scyttools.R --quality_control "/Volumes/Group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF732_190328_G6/outs/filtered_feature_bc_matrix/" "~/lucap1731/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/lucap1731/qc-sce.Rdata" "~/lucap1731/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/lucap1731/dim-red-sce.Rdata" "~/lucap1731/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/lucap1731/clustered-sce.Rdata" "~/lucap1731/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/lucap1731/trajectoried-sce.Rdata" "~/lucap1731/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF731_190328_G18/outs/filtered_feature_bc_matrix/" "~/lucap1452/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/lucap1452/qc-sce.Rdata" "~/lucap1452/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/lucap1452/dim-red-sce.Rdata" "~/lucap1452/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/lucap1452/clustered-sce.Rdata" "~/lucap1452/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/lucap1452/trajectoried-sce.Rdata" "~/lucap1452/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF730_190328_G7/outs/filtered_feature_bc_matrix/" "~/mb44/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/mb44/qc-sce.Rdata" "~/mb44/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/mb44/dim-red-sce.Rdata" "~/mb44/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/mb44/clustered-sce.Rdata" "~/mb44/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/mb44/trajectoried-sce.Rdata" "~/mb44/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group09/CCB/Beshiri/Folders_old/Omics_data/CT35_10x_filtered_gbm/" "~/CT35-2-v1/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/CT35-2-v1/qc-sce.Rdata" "~/CT35-2-v1/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/CT35-2-v1/dim-red-sce.Rdata" "~/CT35-2-v1/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/CT35-2-v1/clustered-sce.Rdata" "~/CT35-2-v1/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/CT35-2-v1/trajectoried-sce.Rdata" "~/CT35-2-v1/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group09/CCB/Beshiri/Folders_old/Omics_data/ct60_filtered_gene_bc_matrices/GRCh38/" "~/ct60/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/ct60/qc-sce.Rdata" "~/ct60/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/ct60/dim-red-sce.Rdata" "~/ct60/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/ct60/clustered-sce.Rdata" "~/ct60/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/ct60/trajectoried-sce.Rdata" "~/ct60/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group09/CCB/Beshiri/Folders_old/Omics_data/mb155_filtered_gene_bc_matrices/GRCh38/" "~/MB155/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/MB155/qc-sce.Rdata" "~/MB155/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/MB155/dim-red-sce.Rdata" "~/MB155/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/MB155/clustered-sce.Rdata" "~/MB155/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "~/MB155/trajectoried-sce.Rdata" "~/MB155/geneset-scored-sce.Rdata"
