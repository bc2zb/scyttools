#!/bin/bash

Rscript scyttools.R --quality_control "/Volumes/Group09/CCB/Beshiri/Folders_old/CT35/'Omics_data/single_cell_RNAseq/lineage-tracing-May-2018/CT35_10x_filtered_gbm" ~/scyttools_test_out/
Rscript scyttools.R --dimensionality_reduction "~/scyttools/scyttools_sce_quality_control.Rdata" ~/scyttools_test_out/
Rscript scyttools.R --cluster_cells "~/scyttools/scyttools_sce_dimensionality_reduction.Rdata" ~/scyttools_test_out/