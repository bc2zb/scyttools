#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
scyttools.R (-h | --help | --version)
scyttools.R --quality_control DIR OUT
scyttools.R --dimensionality_reduction RDS OUT
scyttools.R --cluster_cells RDS OUT
scyttools.R --trajectory_inference RDS OUT
scyttools.R --geneset_scoring RDS OUT
scyttools.R --differential_expression RDS OUT

Description:   This program is a command line interface to running automated single cell RNA sequencing analysis

Options:

--version                     Show the current version.
--quality_control             Perform quality control and remove cells leveraging the Scater and Scran packages
--dimensionality_reduction    Perform dimensionality reduction using GLM PCA
--cluster_cells               Cluster cells using GLM PCA and self-organizing maps
--trajectory_inference        Perform trajectory inference using monocle version 2
--geneset_scoring             Score cells and clusters for geneset activity
--differential_expression     Perform differential expression analysis of clusters

Arguments:

DIR   Directory that contains cellranger gene barcode matrix output
OUT   Provide output file
RDS   RDS file that contains saved SingleCellExperiment object


" -> doc

args <- docopt(doc)

if(args$`--version` == T){ # returns version if version is requested
  cat("\nVersion is 0.1\n")
}else{ # Analysis begins
  # create new sub directory for each instance of cyttools
  RESULTS_DIR <- args$OUT
  dir.create(RESULTS_DIR, showWarnings = F, recursive = T)
  argsFileName <- paste(RESULTS_DIR, "scyttools.args", ".Rdata", sep = "")
  save(args, file = argsFileName)
  
  if(args$`--quality_control` == T){
    
    COMMAND <- paste("Rscript scyttools_quality_control.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--dimensionality_reduction` == T){
    
    COMMAND <- paste("Rscript scyttools_dimensionality_reduction.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster_cells` == T){
    
    COMMAND <- paste("Rscript scyttools_cluster_cells.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--trajectory_inference` == T){
    
    COMMAND <- paste("Rscript scyttools_trajectory_inference.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--geneset_scoring` == T){
    
    COMMAND <- paste("Rscript scyttools_geneset_scoring.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--differential_expression` == T){
    
    COMMAND <- paste("Rscript scyttools_differential_expression.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else{
    cat(paste(c("\nWARNING:","\ncommand not found:", args$`--cluster`), collapse = "\n"), "\n")  
  } 

}
