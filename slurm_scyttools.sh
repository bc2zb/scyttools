#! /bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --partition=ccr
#SBATCH --mem=128g
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --mail-user=capaldobj@nih.gov

module load R

# /data/capaldobj/CS024650_Kelly_Agarwal_2019_07_17/02_PrimaryAnalysisOutput/00_FullCellrangerOutput/Seq1/SCAF604_683/outs/filtered_feature_bc_matrix/
# /data/capaldobj/CS024650_Kelly_Agarwal_2019_07_17/02_PrimaryAnalysisOutput/00_FullCellrangerOutput/Seq1/SCAF605_685/outs/filtered_feature_bc_matrix/

Rscript 
