#' A function to visualise the features distribtuion
#'
#'
#' @param sce A singlecellexperiment object
#' @param RNA_exprs_value A character indicates which expression value for RNA in assayNames is used.
#' @param altExp_name A character indicates which expression matrix is used. by default is none (i.e. RNA).
#' @param altExp_exprs_value A character indicates which expression value in assayNames is used.
#' @param RNA_feature_subset A vector of characters indicates the subset of features of RNA that are used for visualisation
#' @param ADT_feature_subset A vector of characters indicates the subset of features of ADT that are used for visualisation
#' @param cell_subset A vector of characters indicates the subset of cells that are used for visualisation
#' @param cor_threshold Thresholds of correlation.
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param RNA_exprs_pct A numeric indicates the threshold expression percentage of a gene to be considered in correlation analysis
#' @param ADT_exprs_pct A numeric indicates the threshold expression percentage of a gene to be considered in correlation analysis
#' @param RNA_exprs_threshold A numeric indicates the threshold of RNA expression. By default is 0.
#' @param ADT_exprs_threshold A numeric indicates the threshold of ADT expression. By default is 0.
#' @param network_layout layout of the network
#' @param return_igraph indicates whether return the igraph object
#'
#'
#' @importFrom reshape2 melt
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom igraph V E graph_from_data_frame
#' @importFrom graphics legend plot
#' @export

geneADTnetwork <- function(sce,
                           RNA_exprs_value = "logcounts",
                           altExp_name = "ADT",
                           altExp_exprs_value = "logcounts",
                           RNA_feature_subset = NULL,
                           ADT_feature_subset = NULL,
                           cell_subset = NULL,
                           cor_threshold = 0.5,
                           cor_method = c("pearson", "kendall", "spearman"),
                           RNA_exprs_pct = 0.1,
                           ADT_exprs_pct = 0.1,
                           RNA_exprs_threshold = 0,
                           ADT_exprs_threshold = 0,
                           network_layout = NULL,
                           return_igraph = FALSE
) {

  cor_method <- match.arg(cor_method, c("pearson", "kendall", "spearman"))

  if (!is.null(cell_subset)) {
    if (sum(!cell_subset %in% colnames(sce)) != 0) {
      stop("sce does not contain some or all of cell_subset as cell names")
    }
  } else {
    cell_subset <- colnames(sce)
  }




  # RNA exprssion matrix

  if (!RNA_exprs_value %in% SummarizedExperiment::assayNames(sce)) {
    stop("sce does not contain RNA_exprs_value")
  }

  exprsMat1 <- SummarizedExperiment::assay(sce[, cell_subset], RNA_exprs_value)

  if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
    stop("sce does not contain altExp_name as altExpNames")
  }

  if (!altExp_exprs_value %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
    stop("sce does not contain altExp_exprs_value as assayNames for altExp")
  }

  # ADT exprssion matrix
  exprsMat2 <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce[, cell_subset], altExp_name), altExp_exprs_value)


  if (!is.null(RNA_feature_subset)) {
    if (sum(!RNA_feature_subset %in% rownames(exprsMat1)) != 0) {
      stop("sce does not contain some or all of RNA_feature_subset as features")
    }
  }

  if (!is.null(ADT_feature_subset)) {
    if (sum(!ADT_feature_subset %in% rownames(exprsMat2)) != 0) {
      stop("sce does not contain some or all of ADT_feature_subset as features")
    }
  }




  exprsMat1 <- as.matrix(exprsMat1[RNA_feature_subset, , drop = FALSE])
  exprsMat2 <- as.matrix(exprsMat2[ADT_feature_subset, , drop = FALSE])


  # Filter the features that are not expressed at all...

  meanPct1 <- rowMeans(exprsMat1 > RNA_exprs_threshold)
  meanPct2 <- rowMeans(exprsMat2 > ADT_exprs_threshold)

  exprsMat1 <- exprsMat1[meanPct1 > RNA_exprs_pct, , drop = FALSE]
  exprsMat2 <- exprsMat2[meanPct2 > ADT_exprs_pct, , drop = FALSE]


  corMat <- stats::cor(t(exprsMat1), t(exprsMat2), method = cor_method)

  rna_label <- rownames(exprsMat1)
  rna_id <- paste("RNA", rownames(exprsMat1), sep = "_")

  adt_label <- rownames(exprsMat2)
  adt_id <- paste("ADT", adt_label, sep = "_")

  colnames(corMat) <- adt_id
  rownames(corMat) <- rna_id

  df_corMat <- reshape2::melt(corMat)

  df_corMat <- df_corMat[abs(df_corMat$value) > cor_threshold, ]

  g <- igraph::graph_from_data_frame(df_corMat,
                                     directed = FALSE)

  igraph::V(g)$label <- unlist(lapply(strsplit(names(igraph::V(g)), "_"), function(x) paste(x[-1], collapse = "_")))
  igraph::V(g)$class <- unlist(lapply(strsplit(names(igraph::V(g)), "_"), "[[", 1))
  igraph::V(g)$type <- c(TRUE, FALSE)[as.numeric(as.factor(igraph::V(g)$class))]
  igraph::V(g)$shape <- c("circle", "square")[as.numeric(as.factor(igraph::V(g)$class))]
  igraph::V(g)$color <- c("#A0CBE8", "#FFBE7D")[as.numeric(as.factor(igraph::V(g)$class))]
  igraph::V(g)$size <- 10
  igraph::V(g)$label.cex <- 0.4
  igraph::V(g)$label.color <- "black"

  igraph::E(g)$color <-  ifelse(df_corMat$value > 0,
                                "#E15759", "#4E79A7")
  igraph::E(g)$weights <- abs(df_corMat$value) * 10


  graphics::plot(g, layout = network_layout)

  graphics::legend('topleft',
                   legend = c(levels(as.factor(igraph::V(g)$class)), "positive", "negative"),
                   col = c("#A0CBE8", "#FFBE7D", "#E15759", "#4E79A7"),
                   pch = c(16, 15, 95, 95))

  return(g)

}
