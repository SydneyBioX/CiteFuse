#' importanceADT
#'
#' A function to calculate the importance score of ADT
#'
#' @param sce A singlecellexperiment object
#' @param altExp_name A character indicates which expression
#' matrix is used. by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value
#' in assayNames is used.
#' @param method A character indicates the method of ADT importance calculation,
#' either randomForest or PCA
#' @param group A vector indicates the grouping of the data (for random forest)
#' @param subsample Whether perform subsampling (for random forest)
#' @param times A numeric indicates the times of subsampling is performed
#' (for random forest)
#' @param prop A numeric indicates the proportion of cells are subsampled
#' from the whole data (for random forest)
#' @param k_pca Number of principal component will be used to
#' calculate the loading scores (for PCA)
#' @param remove_first_PC A logical input indicates whether
#' the first component will be removed from calculation (for PCA).
#' @param ... other arguments to `randomForest()` or `prcomp()` function
#'
#' @return A SingleCellExperiment object
#'
#' @examples
#' data("sce_control_subset", package = "CiteFuse")
#' sce_control_subset <- importanceADT(sce_control_subset,
#' group = sce_control_subset$SNF_W_louvain,
#' subsample = TRUE)
#'
#' @importFrom randomForest randomForest
#' @importFrom stats prcomp
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom S4Vectors metadata
#'
#' @details
#' For random forest, the importance scores are based on features importance.
#' For PCA, it implements the method proposed in Levin et al
#' (based on the loading of features).
#'
#' @references
#' Levine, J.H., Simonds, E.F., Bendall, S.C., Davis, K.L.,
#' El-ad, D.A., Tadmor, M.D., Litvin, O., Fienberg, H.G., Jager, A.,
#' Zunder, E.R. and Finck, R., 2015.
#' Data-driven phenotypic dissection of AML reveals progenitor-like cells that
#' correlate with prognosis. Cell, 162(1), pp.184-197.
#'
#' @export

importanceADT <- function(sce,
                          altExp_name = "ADT",
                          exprs_value = "logcounts",
                          method = c("randomForest", "PCA"),
                          group = NULL,
                          subsample = TRUE,
                          times = 10,
                          prop = 0.8,
                          k_pca = 5,
                          remove_first_PC = TRUE,
                          ...) {

  if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
    stop("sce does not contain altExp_name as altExpNames")
  }

  if (!exprs_value %in% assayNames(altExp(sce, altExp_name))) {
    stop("sce does not contain exprs_value as assayNames for altExp")
  }

  method <- match.arg(method, c("randomForest", "PCA"))

  if (method == "randomForest" & is.null(group)) {
    stop("To run randomForest ADT importance calculation,
         please provide group infomation.")
  }

  if (k_pca < 2 | k_pca <= 2 & remove_first_PC) {
    stop("Please set a larger k_pca")
  }


  exprsMat <- assay(SingleCellExperiment::altExp(sce, altExp_name),
                    exprs_value)

  if (method == "randomForest") {
    if (subsample) {
      num_sub <- round(ncol(exprsMat) * prop)
      rf <- lapply(seq_len(times), function(x) {
        idx <- sample(ncol(exprsMat), num_sub)
        randomForest::randomForest(t(as.matrix(exprsMat[, idx])),
                                   as.factor(droplevels(group)[idx]), ...)
      })

      importance <- do.call(cbind, lapply(rf, function(x) x$importance))
      colnames(importance) <- seq_len(times)
    } else {
      rf <- randomForest::randomForest(t(as.matrix(exprsMat)),
                                       as.factor(droplevels(group)), ...)
      importance <- rf$importance
    }

    S4Vectors::metadata(sce)[["importanceADT_matrix"]] <- as.matrix(importance)
    S4Vectors::metadata(sce)[["importanceADT"]] <- Matrix::rowMeans(as.matrix(importance))
  }
  if (method == "PCA") {
    adt_pca <- stats::prcomp(t(exprsMat), scale = TRUE, center = TRUE)

    if (remove_first_PC) {
      use_k <- seq(2, k_pca)
    } else {
      use_k <- seq_len(k_pca)
    }
    pca_eigenvalue <- adt_pca$sdev[use_k]^2
    # eigenvalue
    ev_mat <- matrix(1, ncol = 1, nrow = nrow(adt_pca$rotation)) %*%
      matrix(pca_eigenvalue, ncol = length(use_k))

    NRS <- Matrix::rowSums(abs(adt_pca$rotation[, use_k]) * ev_mat)
    S4Vectors::metadata(sce)[["importanceADT_matrix"]] <- as.matrix(NRS)
    S4Vectors::metadata(sce)[["importanceADT"]] <- rowMeans(as.matrix(NRS))
  }


  return(sce)
}


#' visImportance
#'
#' A function to visualise the features distribtuion
#'
#'
#' @param sce A singlecellexperiment object
#' @param plot A string indicates the type of the plot
#'  (either boxplot or heatmap)
#' @param altExp_name A character indicates which expression matrix
#' is used. by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value
#' in assayNames is used.
#'
#' @return A plot (either ggplot or pheatmap) to visualise
#' the ADT importance results
#'
#' @examples
#' data("sce_control_subset", package = "CiteFuse")
#' sce_control_subset <- importanceADT(sce_control_subset,
#' group = sce_control_subset$SNF_W_louvain,
#' subsample = TRUE)
#' visImportance(sce_control_subset, plot = "boxplot")
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
    scores <- scores[order(Matrix::rowMeans(as.matrix(scores))), ,
                     drop = FALSE]



    df_toPlot <- reshape2::melt(scores)
    colnames(df_toPlot) <- c("features", "iter", "importance")

    g <- ggplot(df_toPlot, aes(x = df_toPlot$features,
                               y = df_toPlot$importance,
                               color = df_toPlot$features)) +
      geom_boxplot(outlier.size = 1, outlier.stroke = 0.3,
                   outlier.alpha = 0.8,
                   width = 0.3) +
      scale_colour_viridis_d(direction = -1, end = 0.95) +
      coord_flip() +
      theme_bw() +
      ylab("Importance Scores") +
      xlab("") +
      labs(col = "Features") +
      theme(legend.position = "none")
    print(g)

  }

  if (plot == "heatmap") {


    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
      stop("sce does not contain altExp_name as altExpNames")
    }

    if (!exprs_value %in% assayNames(altExp(sce, altExp_name))) {
      stop("sce does not contain exprs_value as assayNames for altExp")
    }


    corMat_adt <- cor(t(as.matrix(assay(altExp(sce, altExp_name),
                                        exprs_value))))
    corMat_adt[is.na(corMat_adt)] <- 0


    anno_col <- data.frame(importance = rowMeans(scores))
    rownames(anno_col) <- colnames(corMat_adt)
    pheatmap::pheatmap(corMat_adt,
                       clustering_method = "ward.D2",
                       annotation_col = anno_col,
                       breaks = seq(-0.8, 0.8, 1.6/100))
  }



}






