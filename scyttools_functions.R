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

listCols<-function(m){
  #converts a sparse Matrix into a list of its columns
  expr <- split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
  expr
}

score_signature <- function(m, sig_weights){
  expr <- split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
  anno <- split(rownames(m)[summary(m)$i], findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
  scores <- lapply(seq_along(expr), function(x){
    cell_expr <- data.frame(expr = expr[[x]],
                            ensgene = anno[[x]]) %>% 
      inner_join(sig_weights,
                 by = "ensgene")
    cell_score <- sum(cell_expr$expr*cell_expr$V1)
    return(cell_score)
  })
  names(scores) <- colnames(m)
  return(scores)
}


# baseline3<-function(m,f){
#   vapply(listCols(m), f, FUN.VALUE=0.0)
# }
# 
# f<-median
# 
# baseline3(logcounts(sce_glm_pca),f)
# 
# 
# expr_values <- listCols(play_with_me)  
