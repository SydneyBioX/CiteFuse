#' visualiseExprs
#'
#' A function to visualise the features distribtuion
#'
#'
#' @param sce A singlecellexperiment object
#' @param plot Type of plot, includes boxplot, violin, jitter, density, and pairwise. By default is boxplot
#' @param altExp_name A character indicates which expression matrix is used. by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value in assayNames is used.
#' @param feature_subset A vector of characters indicates the subset of features that are used for visualisation
#' @param cell_subset A vector of characters indicates the subset of cells that are used for visualisation
#' @param n A numeric indicates the top expressed features to show.
#' @param threshold Thresholds of high expresion for features (only is used for pairwise plot).
#'
#' @return A ggplot to visualise te features distribution
#'
#' @importFrom Matrix rowMeans
#' @importFrom reshape2 melt
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom  ggridges geom_density_ridges2
#' @importFrom gridExtra grid.arrange
#' @importFrom utils combn
#' @import ggplot2
#' @export

visualiseExprs <- function(sce,
                           plot = c("boxplot", "violin", "jitter", "density", "pairwise"),
                           altExp_name = c("none"),
                           exprs_value = "logcounts",
                           feature_subset = NULL,
                           cell_subset = NULL,
                           n = NULL,
                           threshold = NULL
                           ) {

  plot <- match.arg(plot, c("boxplot", "violin", "jitter", "density", "pairwise"))

  if (!is.null(cell_subset)) {
    if (sum(!cell_subset %in% colnames(sce)) != 0) {
      stop("sce does not contain some or all of cell_subset as cell names")
    }
  } else {
    cell_subset <- colnames(sce)
  }


  if (altExp_name != "none") {
    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
      stop("sce does not contain altExp_name as altExpNames")
    }

    if (!exprs_value %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
      stop("sce does not contain exprs_value as assayNames for altExp")
    }

    exprsMat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce[, cell_subset], altExp_name), exprs_value)

  } else {

    # if altExp_name is "none", then the assay in SingleCellExperiment is extracted (RNA in most of the cases)

    exprsMat <- SummarizedExperiment::assay(sce[, cell_subset], exprs_value)
  }





  if (is.null(feature_subset)) {
    if (is.null(n)) {
      n <- min(nrow(exprsMat), 10)
    }
    rowmeans <- apply(exprsMat, 1, median)
    rowmeans <- sort(rowmeans, decreasing = TRUE)
    feature_subset <- rev(names(rowmeans)[seq_len(n)])
  } else {
    if (sum(!feature_subset %in% rownames(exprsMat)) != 0) {
      stop("sce does not contain some or all of feature_subset as features")
    }
  }

  exprsMat <- as.matrix(exprsMat[feature_subset, , drop = FALSE])

  if (plot == "pairwise") {
    combination <- utils::combn(feature_subset, 2)

    ggList <- apply(combination, 2, function(x) {
      suppressMessages(scatterSingle(exprsMat[x, ], threshold = threshold))
    })



    do.call(gridExtra::grid.arrange, c(ggList, ncol = min(length(ggList), 2)))
  } else {
    df_toPlot <- reshape2::melt(exprsMat)
    colnames(df_toPlot) <- c("features", "cells", "value")


    if (plot == "boxplot") {
      g <- ggplot(df_toPlot, aes(x = df_toPlot$features,
                                 y = df_toPlot$value,
                                 color = df_toPlot$features)) +
        geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                     width = 0.3) +
        scale_colour_viridis_d(direction = -1, end = 0.95) +
        coord_flip() +
        theme_bw() +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }

    if (plot == "violin") {
      g <- ggplot(df_toPlot, aes(x = df_toPlot$features, y = df_toPlot$value)) +
        geom_violin(aes(fill = df_toPlot$features)) +
        geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                     width = 0.05) +
        scale_fill_viridis_d(direction = -1, end = 0.95) +
        coord_flip() +
        theme_bw() +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }

    if (plot == "jitter")  {
      g <- ggplot(df_toPlot, aes(x = df_toPlot$features, y = df_toPlot$value)) +
        geom_jitter(aes(color = df_toPlot$features), size = 0.5, alpha = 0.5) +
        scale_color_viridis_d(direction = -1, end = 0.95) +
        coord_flip() +
        theme_bw() +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }

    if (plot == "density")  {
      g <- ggplot(df_toPlot, aes(x = df_toPlot$value, y = df_toPlot$features)) +
        ggridges::geom_density_ridges2(aes(fill = df_toPlot$features), alpha = 0.5) +
        scale_fill_viridis_d(direction = -1, end = 0.95) +
        theme_bw() +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }
    return(g)
  }

}







# A function to generate scatter plot pairwise feature expression


scatterSingle <- function(exprsMat, group = NULL, threshold = NULL) {


  df <- data.frame(t(as.matrix(exprsMat)))

  colnames(df) <- lapply(strsplit(rownames(exprsMat), "\\."), "[[", 1)


  if (is.null(group)) {

    if (is.null(threshold)) {
      threshold <- apply(exprsMat, 1, fitMixtures)
    }



    exprsMat_pass <- exprsMat > threshold

    nn <- paste(c(rownames(exprsMat), ""), collapse = "-")
    pp <- paste(c(rownames(exprsMat), ""), collapse = "+")
    np <- paste(rownames(exprsMat)[1], "-", rownames(exprsMat)[2], "+", sep = "")
    pn <- paste(rownames(exprsMat)[1], "+", rownames(exprsMat)[2], "-", sep = "")
    group <- rep(nn , ncol(exprsMat))
    group[exprsMat_pass[1, ] & exprsMat_pass[2, ]] <- pp
    group[exprsMat_pass[1, ] & !exprsMat_pass[2, ]] <- pn
    group[!exprsMat_pass[1, ] & exprsMat_pass[2, ]] <- np

    group <- factor(group, levels = c(nn, np, pn, pp))

    point_col <- c("grey50", "#377EB8", "#E41A1C", "#984EA3")

    names(point_col) <- levels(group)

    col_layer <- scale_color_manual(values = point_col)
  } else{
    col_layer <- scale_color_hue()
  }



  pmain <- ggplot(df, aes(x = df[, 1], y = df[, 2], color = group)) +
    geom_point(alpha = 0.5) +
    # geom_density2d(bins = 10) +
    # stat_density2d(aes(alpha = ..level.., fill = ..level..), geom = "polygon", n = 100, bins = 10) +
    geom_hline(yintercept = threshold[2], col = "red", linetype = 2, size = 1) +
    geom_vline(xintercept = threshold[1], col = "red", linetype = 2, size = 1) +
    col_layer +
    scale_x_continuous(position = "top") +
    scale_y_continuous(position = "right") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2])




  xdens <- cowplot::axis_canvas(pmain, axis = "x") +
    geom_density(data = df, aes(x = df[, 1], fill = exprsMat_pass[1, ], alpha = 0.5)) +
    scale_fill_manual(values = c("grey50", "#E41A1C")) +
    scale_y_reverse() +
    NULL

  ydens <- cowplot::axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
    geom_density(data = df, aes(x = df[, 2], fill = exprsMat_pass[2, ], alpha = 0.5)) +
    coord_flip() +
    scale_y_reverse() +
    scale_fill_manual(values = c("grey50", "#377EB8")) +
    NULL

  p1 <- cowplot::insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "bottom")

  p2 <- cowplot::insert_yaxis_grob(p1, ydens,
                                   grid::unit(.2, "null"), position = "left")

  # ggdraw(p2)
  return(p2)
}


fitMixtures <- function(vec) {
  km <- stats::kmeans(vec, centers = 2, nstart = 1000)
  mixmdl <- try(mixtools::normalmixEM(vec,
                                      fast = TRUE, maxrestarts = 1000,
                                      k = 2, maxit = 10000,
                                      mu = km$centers[,1],
                                      ECM = TRUE, verb = FALSE),
                silent = TRUE)



  if (class(mixmdl) == "try-error") {
    threshold <- min(max(vec[km$cluster == 1]), max(vec[km$cluster == 2]))
  } else {
    threshold <- getThreshold(mixmdl)
  }
  return(threshold)
}




