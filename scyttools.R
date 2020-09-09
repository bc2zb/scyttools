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
scyttools.R --export_results RDS OUT
scyttools.R integrate_batchelor OUT RDS RDS...


Description:   This program is a command line interface to running automated single cell RNA sequencing analysis

Options:

--version                     Show the current version.
--quality_control             Perform quality control and remove cells leveraging the Scater and Scran packages
--dimensionality_reduction    Perform dimensionality reduction using GLM PCA
--cluster_cells               Cluster cells using GLM PCA and self-organizing maps
--trajectory_inference        Perform trajectory inference using monocle version 3
--geneset_scoring             Score cells and clusters for geneset activity
--export_results              Perform differential expression analysis of clusters and export
--integrate_batchelor         Integrates multiple runs together

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
  RESULTS_DIR <- dirname(args$OUT)
  if(RESULTS_DIR != "."){
    dir.create(RESULTS_DIR, showWarnings = F, recursive = T)
  }
  
  if(args$`--quality_control` == T){
    
    COMMAND <- paste("Rscript scyttools_quality_control.R",
                     paste("'", args$DIR, "'", sep = ""),
                     paste("'", args$OUT, "'", sep = ""))
    
    cat(COMMAND, "\n")
    
    system(command = COMMAND)
  }else if(args$`--dimensionality_reduction` == T){
    
    COMMAND <- paste("Rscript scyttools_dimensionality_reduction.R",
                     paste("'", args$RDS, "'", sep = ""),
                     paste("'", args$OUT, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster_cells` == T){
    
    COMMAND <- paste("Rscript scyttools_cluster_cells.R",
                     paste("'", args$RDS, "'", sep = ""),
                     paste("'", args$OUT, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--trajectory_inference` == T){
    
    COMMAND <- paste("Rscript scyttools_trajectory_inference.R",
                     paste("'", args$RDS, "'", sep = ""),
                     paste("'", args$OUT, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--geneset_scoring` == T){
    
    COMMAND <- paste("Rscript scyttools_geneset_scoring.R",
                     paste("'", args$RDS, "'", sep = ""),
                     paste("'", args$OUT, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`integrate_batchelor` == T){
    
    COMMAND <- paste("Rscript scyttools_integrate.R",
                     paste("'", args$OUT, "'", sep = ""),
                     paste("'", 
                           paste(args$RDS, collapse = "' '"),
                           "'",
                           sep = "")
    )
    system(command = COMMAND)
    
  }else if(args$`--export_results` == T){
    
    COMMAND <- paste("Rscript scyttools_export_results.R",
                     paste("'", args$RDS, "'", sep = ""),
                     paste("'", args$OUT, "'", sep = ""))
    system(command = COMMAND)
  }else{
    cat(paste(c("\nWARNING:","\ncommand not found:"), collapse = "\n"), "\n")  
  } 

}
