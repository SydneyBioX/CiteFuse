#' A function to calculate the importance score of ADT
#' 
#' 
#' @param sce A singlecellexperiment object
#' @param altExp_name A character indicates which expression matrix is used. by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value in assayNames is used.
#' @param group A vector indicates the grouping of the data
#' @param subsample Whether perform subsampling
#' @param times A numeric indicates the times of subsampling is performed
#' @param prop A numeric indicates the proportion of cells are subsampled from the whole data
#' @param ... other arguments to `randomForest()` function
#' 
#' @importFrom randomForest randomForest
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom S4Vectors metadata
#' 
#' @export

importanceADT <- function(sce, 
                          altExp_name = "ADT", 
                          exprs_value = "logcounts",
                          group = NULL,
                          subsample = TRUE,
                          times = 10,
                          prop = 0.8,
                          ...) {
  
  if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
    stop("sce does not contain altExp_name as altExpNames")
  }
  
  if (!exprs_value %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
    stop("sce does not contain exprs_value as assayNames for altExp")
  }
  
  
  
  exprsMat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), exprs_value)
  
  if (subsample) {
    num_sub <- round(ncol(exprsMat) * prop)
    rf <- lapply(seq_len(times), function(x) {
      idx <- sample(ncol(exprsMat), num_sub)
      randomForest::randomForest(t(as.matrix(exprsMat[, idx])), as.factor(group[idx]), ...)
    })
    
    importance <- do.call(cbind, lapply(rf, function(x) x$importance))
    colnames(importance) <- seq_len(times)
  } else {
    rf <- randomForest::randomForest(t(as.matrix(exprsMat)), as.factor(group), ...)
    importance <- rf$importance
  }
  
  S4Vectors::metadata(sce)[["importanceADT_matrix"]] <- as.matrix(importance)
  S4Vectors::metadata(sce)[["importanceADT"]] <- Matrix::rowMeans(as.matrix(importance))
  
  return(sce)
}


#' A function to visualise the features distribtuion
#' 
#' 
#' @param sce A singlecellexperiment object
#' @param plot A string indicates the type of the plot (either boxplot or heatmap)
#' @param altExp_name A character indicates which expression matrix is used. by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value in assayNames is used.
#' 
#' @importFrom Matrix rowMeans
#' @importFrom reshape2 melt
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom S4Vectors metadata
#' @importFrom pheatmap pheatmap
#' @import ggplot2
#' @export

visImportance <- function(sce, 
                          plot = c("boxplot", "heatmap"),
                          altExp_name = "ADT", 
                          exprs_value = "logcounts") {
  
  
  if (!"importanceADT" %in% names(S4Vectors::metadata(sce))) {
    stop("Please perform importanceADT() first")
  }
  
  
  
  plot <- match.arg(plot, c("boxplot", "heatmap"))
  
  scores <- metadata(sce)[["importanceADT_matrix"]]
  
  
  if (plot == "boxplot") {
    scores <- scores[order(Matrix::rowMeans(as.matrix(scores))), , drop = FALSE]
    
    
    
    df_toPlot <- reshape2::melt(scores)
    colnames(df_toPlot) <- c("features", "iter", "importance")
    
    g <- ggplot(df_toPlot, aes(x = features, y = importance, color = features)) +
      geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                   width = 0.3) +
      scale_colour_viridis_d(direction = -1, end = 0.95) +
      coord_flip() +
      theme_bw() +
      ylab("Importance Scores") +
      xlab("") +
      theme(legend.position = "none")
    print(g)
    
  }
  
  if (plot == "heatmap") {
    
    
    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
      stop("sce does not contain altExp_name as altExpNames")
    }
    
    if (!exprs_value %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
      stop("sce does not contain exprs_value as assayNames for altExp")
    }
    
    
    corMat_adt <- cor(t(as.matrix(SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), exprs_value))))
    
    
    anno_col <- data.frame(importance = rowMeans(scores))
    rownames(anno_col) <- colnames(corMat_adt)
    pheatmap::pheatmap(corMat_adt,
                       clustering_method = "ward.D2",
                       annotation_col = anno_col,
                       breaks = seq(-0.8, 0.8, 1.6/100))
  }

  
  
}






