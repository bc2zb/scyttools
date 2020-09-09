library(DropletUtils)
library(scater)
library(scran)
library(annotables)
library(tidyverse)
library(glmpca)
library(kohonen)
library(ggraph)
library(tidygraph)
library(igraph)
library(Seurat)
library(edgeR)
library(monocle3)
library(scDblFinder)
library(scry)
library(batchelor)

set.seed(8675309)

cor_dist <- function(som_codes, method){
  corr_mat <- cor(t(som_codes),
                  method = method)
  return(1 - corr_mat)
}

tidy_sparse_matrix <- function (x) 
{
  s <- Matrix::summary(x)
  row <- s$i
  if (!is.null(rownames(x))) {
    row <- rownames(x)[row]
  }
  col <- s$j
  if (!is.null(colnames(x))) {
    col <- colnames(x)[col]
  }
  ret <- data.frame(row = row, column = col, value = s$x, stringsAsFactors = FALSE)
  ret
}
