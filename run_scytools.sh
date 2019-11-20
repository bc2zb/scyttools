#!/bin/bash

#Rscript scyttools.R --quality_control "/Volumes/Group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF732_190328_G6/outs/filtered_feature_bc_matrix/" "~/lucap1731/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/lucap1731/qc-sce.Rdata" "~/lucap1731/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "~/lucap1731/dim-red-sce.Rdata" "~/lucap1731/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "~/lucap1731/clustered-sce.Rdata" "~/lucap1731/trajectoried-sce.Rdata"
#Rscript scyttools.R --geneset_scoring "~/lucap1731/trajectoried-sce.Rdata" "~/lucap1731/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF731_190328_G18/outs/filtered_feature_bc_matrix/" "~/lucap1452/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/lucap1452/qc-sce.Rdata" "~/lucap1452/dim-red-sce.Rdata"
#Rscript scyttools.R --cluster_cells "~/lucap1452/dim-red-sce.Rdata" "~/lucap1452/clustered-sce.Rdata"
#Rscript scyttools.R --trajectory_inference "~/lucap1452/clustered-sce.Rdata" "~/lucap1452/trajectoried-sce.Rdata"
#Rscript scyttools.R --geneset_scoring "~/lucap1452/trajectoried-sce.Rdata" "~/lucap1452/geneset-scored-sce.Rdata"

#Rscript scyttools.R --quality_control "/Volumes/Group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF730_190328_G7/outs/filtered_feature_bc_matrix/" "~/mb44/qc-sce.Rdata"
#Rscript scyttools.R --dimensionality_reduction "~/mb44/qc-sce.Rdata" "~/mb44/dim-red-sce.Rdata"
#Rscript scyttools.R --cluster_cells "~/mb44/dim-red-sce.Rdata" "~/mb44/clustered-sce.Rdata"
#Rscript scyttools.R --trajectory_inference "~/mb44/clustered-sce.Rdata" "~/mb44/trajectoried-sce.Rdata"
#Rscript scyttools.R --geneset_scoring "~/mb44/trajectoried-sce.Rdata" "~/mb44/geneset-scored-sce.Rdata"
