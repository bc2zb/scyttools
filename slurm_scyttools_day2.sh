#! /bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --partition=ccr
#SBATCH --mem=128g
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --mail-user=capaldobj@nih.gov

source ~/.bash_profile
conda activate radian

# /data/capaldobj/CS024650_Kelly_Agarwal_2019_07_17/02_PrimaryAnalysisOutput/00_FullCellrangerOutput/Seq1/SCAF604_683/outs/filtered_feature_bc_matrix/
# /data/capaldobj/CS024650_Kelly_Agarwal_2019_07_17/02_PrimaryAnalysisOutput/00_FullCellrangerOutput/Seq1/SCAF605_685/outs/filtered_feature_bc_matrix/

Rscript scyttools.R --quality_control "/data/capaldobj/cell-tagging-2020-may/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF1385_145_2-CTL-D2/outs/filtered_feature_bc_matrix" "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/qc-sce.Rdata"
Rscript scyttools.R --dimensionality_reduction "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/qc-sce.Rdata" "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/dim-red-sce.Rdata"
Rscript scyttools.R --cluster_cells "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/dim-red-sce.Rdata" "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/clustered-sce.Rdata"
Rscript scyttools.R --trajectory_inference "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/clustered-sce.Rdata" "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/trajectoried-sce.Rdata"
Rscript scyttools.R --geneset_scoring "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/trajectoried-sce.Rdata" "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/geneset-scored-sce.Rdata"
Rscript scyttools.R --export_results "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/trajectoried-sce.Rdata" "/data/capaldobj/cell-tagging-2020-may/lucap_145_2_2_day/"
