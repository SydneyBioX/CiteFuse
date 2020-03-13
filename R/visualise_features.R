#' visualiseExprs
#'
#' A function to visualise the features distribtuion
#'
#'
#' @param sce A singlecellexperiment object
#' @param plot Type of plot, includes boxplot, violin, jitter, density,
#' and pairwise. By default is boxplot
#' @param altExp_name A character indicates which expression matrix is used.
#' by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value
#' in assayNames is used.
#' @param group_by A character indicates how is the expression
#' will be group in the plots (stored in colData).
#' @param facet_by A character indicates how is the expression
#' will be lay out panels in a grid in the plots (stored in colData).
#' @param feature_subset A vector of characters indicates
#' the subset of features that are used for visualisation
#' @param cell_subset A vector of characters indicates
#' the subset of cells that are used for visualisation
#' @param n A numeric indicates the top expressed features to show.
#' @param threshold Thresholds of high expresion for features
#' (only is used for pairwise plot).
#'
#'
#' @return A ggplot to visualise te features distribution
#'
#' @examples
#' data("sce_control_subset", package = "CiteFuse")
#' visualiseExprs(sce_control_subset,
#' plot = "boxplot",
#' group_by = "SNF_W_louvain",
#' feature_subset = c("hg19_CD8A"))
#'
#' visualiseExprs(sce_control_subset,
#' plot = "density",
#' altExp_name = "ADT",
#' group_by = "SNF_W_louvain",
#' feature_subset = c("CD8", "CD4"))
#'
#' @importFrom Matrix rowMeans
#' @importFrom reshape2 melt
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay colData
#' @importFrom  ggridges geom_density_ridges2
#' @importFrom gridExtra grid.arrange
#' @importFrom utils combn
#' @import ggplot2
#' @export

visualiseExprs <- function(sce,
                           plot = c("boxplot", "violin", "jitter",
                                    "density", "pairwise"),
                           altExp_name = c("none"),
                           exprs_value = "logcounts",
                           group_by = NULL,
                           facet_by = NULL,
                           feature_subset = NULL,
                           cell_subset = NULL,
                           n = NULL,
                           threshold = NULL) {

  plot <- match.arg(plot, c("boxplot", "violin", "jitter",
                            "density", "pairwise"))


  if (!is.null(cell_subset)) {
    if (sum(!cell_subset %in% colnames(sce)) != 0) {
      stop("sce does not contain some or all of cell_subset as cell names")
    }
  } else {
    cell_subset <- colnames(sce)
  }


  # if (altExp_name != "none") {
  #   if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
  #     stop("sce does not contain altExp_name as altExpNames")
  #   }
  #
  #   if (!exprs_value %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
  #     stop("sce does not contain exprs_value as assayNames for altExp")
  #   }
  #
  #   exprsMat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce[, cell_subset], altExp_name), exprs_value)
  #
  # } else {
  #
  #   # if altExp_name is "none", then the assay in SingleCellExperiment is extracted (RNA in most of the cases)
  #
  #   exprsMat <- SummarizedExperiment::assay(sce[, cell_subset], exprs_value)
  # }
  #

  exprsMat <- .extract_exprsMat(sce, cell_subset, altExp_name, exprs_value)

  group_by_info <- .extract_group_by(sce, group_by, cell_subset)

  facet_by_info <- .extract_facet_by(sce, facet_by)


  exprsMat <- .extract_filter_features(exprsMat, feature_subset, n)

  exprsMat <- as.matrix(exprsMat)




  if (plot == "pairwise") {



    if (!is.null(threshold)) {
      if (length(threshold) == 1) {
        threshold <- rep(threshold, length(feature_subset))
      }

      if (length(threshold) != length(feature_subset)) {
        stop("length of threshold is not equal to the length of the features subset.")
      }

      names(threshold) <- feature_subset
    }

    combination <- utils::combn(feature_subset, 2)

    ggList <- apply(combination, 2, function(x) {
      suppressMessages(scatterSingle(exprsMat[x, ], threshold = threshold[x]))
    })



    do.call(gridExtra::grid.arrange, c(ggList, ncol = min(length(ggList), 2)))
  } else {

    df_toPlot <- reshape2::melt(exprsMat)
    colnames(df_toPlot) <- c("features", "cells", "value")
    df_toPlot$group <- group_by_info[as.character(df_toPlot$cells)]
    df_toPlot$facet <- facet_by_info[as.character(df_toPlot$cells)]

    if (is.null(facet_by)) {
      g_facet <- facet_grid(~df_toPlot$features)
    } else {
      g_facet <- facet_grid(df_toPlot$features~df_toPlot$facet)
    }

    if (plot == "boxplot") {


      if (!is.null(group_by)) {

        g <- ggplot(df_toPlot, aes(x = group,
                                   y = value,
                                   fill = group)) +
          geom_boxplot(outlier.size = 1, outlier.stroke = 0.3,
                       outlier.alpha = 0.8, width = 0.3) +
          scale_fill_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
          theme_bw() +
          ylab(exprs_value) +
          g_facet +
          xlab("") +
          labs(fill = group_by)

      } else {
        g <- ggplot(df_toPlot, aes(x = features,
                                   y = value,
                                   fill = features)) +
          geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                       width = 0.3) +
          scale_fill_viridis_d(direction = -1, end = 0.95) +
          theme_bw() +
          ylab(exprs_value) +
          xlab("") +
          theme(legend.position = "none")
      }

    }

    if (plot == "violin") {
      if (!is.null(group_by)) {

        g <- ggplot(df_toPlot, aes(x = group,
                                   y = value)) +
          geom_violin(aes(fill = group)) +
          geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                       width = 0.05) +
          scale_fill_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
          theme_bw() +
          ylab(exprs_value) +
          g_facet +
          xlab("") +
          labs(fill = group_by)

      } else {
        g <- ggplot(df_toPlot, aes(x = features,
                                   y = value)) +
          geom_violin(aes(fill = features)) +
          geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                       width = 0.05) +
          scale_fill_viridis_d(direction = -1, end = 0.95) +
          theme_bw() +
          ylab(exprs_value) +
          xlab("") +
          theme(legend.position = "none")
      }
    }

    if (plot == "jitter")  {
      if (!is.null(group_by)) {

        g <- ggplot(df_toPlot, aes(x = group,
                                   y = value)) +
          geom_jitter(aes(color = group), size = 0.5, alpha = 0.5) +
          scale_color_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
          theme_bw() +
          ylab(exprs_value) +
          g_facet +
          xlab("") +
          labs(color = group_by)

      } else {
        g <- ggplot(df_toPlot, aes(x = features, y = value)) +
          geom_jitter(aes(color = features), size = 0.5, alpha = 0.5) +
          scale_color_viridis_d(direction = -1, end = 0.95) +
          coord_flip() +
          theme_bw() +
          ylab(exprs_value) +
          xlab("") +
          theme(legend.position = "none")
      }
    }

    if (plot == "density")  {
      if (!is.null(group_by)) {

        g <- ggplot(df_toPlot, aes(y = group,
                                   x = value)) +
          ggridges::geom_density_ridges2(aes(fill = group),
                                         alpha = 0.5) +
          scale_fill_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
          theme_bw() +
          ylab(group_by) +
          g_facet +
          xlab(exprs_value) +
          labs(fill = group_by)

      } else {
        g <- ggplot(df_toPlot, aes(x = value,
                                   y = features)) +
          ggridges::geom_density_ridges2(aes(fill = features),
                                         alpha = 0.5) +
          scale_fill_viridis_d(direction = -1, end = 0.95) +
          theme_bw() +
          ylab(exprs_value) +
          xlab("") +
          theme(legend.position = "none")
      }
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

    threshold_matrix <- t(sapply(threshold, rep, times = ncol(exprsMat)))

    exprsMat_pass <- exprsMat > threshold_matrix

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



#' visualiseExprsList
#'
#' A function to visualise the features distribtuion for
#' a list of SingleCellExperiment
#'
#'
#' @param sce_list A list of SingleCellExperiment object
#' @param plot Type of plot, includes boxplot, violin, jitter, density,
#' and pairwise. By default is boxplot
#' @param altExp_name A character indicates which expression matrix is used.
#' by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value
#' in assayNames is used.
#' @param group_by A character indicates how is the expression
#' will be group in the plots (stored in colData).
#' @param feature_subset A vector of characters indicates
#' the subset of features that are used for visualisation
#' @param cell_subset A vector of characters indicates
#' the subset of cells that are used for visualisation
#' @param n A numeric indicates the top expressed features to show.
#'
#' @return A ggplot to visualise te features distribution
#'
#' @examples
#' data("sce_control_subset", package = "CiteFuse")
#' data("sce_ctcl_subset", package = "CiteFuse")
#' visualiseExprsList(sce_list = list(control = sce_control_subset,
#' ctcl = sce_ctcl_subset),
#' plot = "boxplot",
#' altExp_name = "none",
#' exprs_value = "logcounts",
#' feature_subset = c("hg19_CD8A"),
#' group_by = c("SNF_W_louvain", "SNF_W_louvain"))
#'
#'
#' @importFrom Matrix rowMeans
#' @importFrom reshape2 melt
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay colData
#' @importFrom  ggridges geom_density_ridges2
#' @importFrom gridExtra grid.arrange
#' @importFrom utils combn
#' @importFrom  methods is
#' @import ggplot2
#' @export

visualiseExprsList <- function(sce_list,
                               plot = c("boxplot", "violin", "jitter",
                                        "density"),
                               altExp_name = "none",
                               exprs_value = "logcounts",
                               group_by = NULL,
                               feature_subset = NULL,
                               cell_subset = NULL,
                               n = NULL) {

  plot <- match.arg(plot, c("boxplot", "violin",
                            "jitter", "density"))

  dataset_name <- names(sce_list)

  if (is.null(dataset_name)) {
    dataset_name <- paste("dataset", seq_len(sce_list), sep = "_")
  }

  if (!is.null(cell_subset)) {

    if (!"list" %in% methods::is(cell_subset) |
        length(cell_subset) != length(sce_list)) {
      stop("cell_subset needs to be a list that
           equal to the length of sce_list")
    }

    for (i in seq_len(length(sce_list))) {
      if (sum(!cell_subset[[i]] %in% colnames(sce_list[[i]])) != 0) {
        stop("sce does not contain some or all of cell_subset as cell names")
      }
    }
  } else{
    cell_subset <- lapply(sce_list, colnames)
  }

  exprsMatList <- lapply(seq_len(length(sce_list)), function(i) {
    exprs <- .extract_exprsMat(sce_list[[i]],
                               cell_subset[[i]],
                               altExp_name, exprs_value)
  })

  names(exprsMatList) <- dataset_name

  if (length(group_by) == 1) {
    group_by <- rep(group_by, length(sce_list))
  }

  group_by_info_list <- lapply(seq_len(length(sce_list)), function(i) {
    group_info <- .extract_group_by(sce_list[[i]],
                                    group_by = group_by[i],
                                    cell_subset = cell_subset[[i]],
                                    factor = FALSE)

    if (!is.null(group_info)) {
      names(group_info) <- paste(colnames(exprsMatList[[i]]),
                                 dataset_name[i], sep = "_")
    }

    group_info
  })

  group_by_info_list <- unlist(group_by_info_list)

  #facet_by_info <- rep(dataset_name, sapply(sce_list, ncol))

  exprsMatList <- lapply(exprsMatList, .extract_filter_features,
                         feature_subset = feature_subset, n = n)



  df_toPlot <- reshape2::melt(exprsMatList)
  colnames(df_toPlot) <- c("features", "cells", "value", "dataset")
  df_toPlot$id <- paste(df_toPlot$cells, df_toPlot$dataset, sep = "_")
  df_toPlot$facet <- df_toPlot$dataset

  if (!is.null(group_by)) {
    df_toPlot$group <- factor(group_by_info_list[df_toPlot$id])
    g_facet <- facet_grid(df_toPlot$features~df_toPlot$facet,
                          scales = "free_x")
  } else {
    g_facet <- facet_grid(~df_toPlot$facet,
                          scales = "free_x")
  }



  if (plot == "boxplot") {

    if (!is.null(group_by)) {

      g <- ggplot(df_toPlot, aes(x = group,
                                 y = value,
                                 fill = group)) +
        geom_boxplot(outlier.size = 1, outlier.stroke = 0.3,
                     outlier.alpha = 0.8, width = 0.3) +
        scale_fill_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
        theme_bw() +
        ylab(exprs_value) +
        g_facet +
        xlab("") +
        labs(fill = group_by) +
        theme(legend.position = "none")

    } else {
      g <- ggplot(df_toPlot, aes(x = features,
                                 y = value,
                                 fill = features)) +
        geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                     width = 0.3) +
        scale_fill_viridis_d(direction = -1, end = 0.95) +
        theme_bw() +
        ylab(exprs_value) +
        g_facet +
        xlab("") +
        theme(legend.position = "none")
    }

  }

  if (plot == "violin") {
    if (!is.null(group_by)) {

      g <- ggplot(df_toPlot, aes(x = group,
                                 y = value)) +
        geom_violin(aes(fill = group)) +
        geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                     width = 0.05) +
        scale_fill_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
        theme_bw() +
        ylab(exprs_value) +
        g_facet +
        xlab("") +
        labs(fill = group_by) +
        theme(legend.position = "none")

    } else {
      g <- ggplot(df_toPlot, aes(x = features,
                                 y = value)) +
        geom_violin(aes(fill = features)) +
        geom_boxplot(outlier.size = 1, outlier.stroke = 0.3, outlier.alpha = 0.8,
                     width = 0.05) +
        scale_fill_viridis_d(direction = -1, end = 0.95) +
        theme_bw() +
        g_facet +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }
  }

  if (plot == "jitter")  {
    if (!is.null(group_by)) {

      g <- ggplot(df_toPlot, aes(x = group,
                                 y = value)) +
        geom_jitter(aes(color = group), size = 0.5, alpha = 0.5) +
        scale_color_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
        theme_bw() +
        ylab(exprs_value) +
        g_facet +
        xlab("") +
        labs(color = group_by) +
        theme(legend.position = "none")

    } else {
      g <- ggplot(df_toPlot, aes(x = features, y = value)) +
        geom_jitter(aes(color = features), size = 0.5, alpha = 0.5) +
        scale_color_viridis_d(direction = -1, end = 0.95) +
        coord_flip() +
        theme_bw() +
        g_facet +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }
  }

  if (plot == "density")  {
    if (!is.null(group_by)) {

      g <- ggplot(df_toPlot, aes(y = group,
                                 x = value)) +
        ggridges::geom_density_ridges2(aes(fill = group),
                                       alpha = 0.5) +
        scale_fill_manual(values = cite_colorPal(nlevels(df_toPlot$group))) +
        theme_bw() +
        ylab(group_by) +
        facet_grid(df_toPlot$features~df_toPlot$facet,
                   scales = "free_y") +
        xlab(exprs_value) +
        labs(fill = group_by) +
        theme(legend.position = "none")

    } else {
      g <- ggplot(df_toPlot, aes(x = value,
                                 y = features)) +
        ggridges::geom_density_ridges2(aes(fill = features),
                                       alpha = 0.5) +
        scale_fill_viridis_d(direction = -1, end = 0.95) +
        theme_bw() +
        g_facet +
        ylab(exprs_value) +
        xlab("") +
        theme(legend.position = "none")
    }
  }
  return(g)


}

#' @importFrom SingleCellExperiment altExpNames
#' @importFrom SummarizedExperiment assayNames assay

.extract_exprsMat <- function(sce, cell_subset, altExp_name, exprs_value) {

  if (!altExp_name %in% c("RNA", "none")) {
    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
      stop("sce does not contain altExp_name as altExpNames")
    }

    if (!exprs_value %in% SummarizedExperiment::assayNames(altExp(sce, altExp_name))) {
      stop("sce does not contain exprs_value as assayNames for altExp")
    }

    exprsMat <- SummarizedExperiment::assay(altExp(sce[, cell_subset],
                                                   altExp_name), exprs_value)

  } else {

    # if altExp_name is "none", then the assay in
    # SingleCellExperiment is extracted (RNA in most of the cases)

    exprsMat <- SummarizedExperiment::assay(sce[, cell_subset], exprs_value)
  }

  return(exprsMat)

}



#' @importFrom SummarizedExperiment colData

.extract_group_by <- function(sce, group_by, cell_subset, factor = TRUE) {
  sce <- sce[, cell_subset]
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(SummarizedExperiment::colData(sce))) {
      stop("group_by is not a column name of sce's colData")
    }
    group_by_info <- as.factor(colData(sce)[, group_by])
    if (!factor) {
      group_by_info <- as.character(group_by_info)
    }
    names(group_by_info) <- colnames(sce)
  } else {
    group_by_info <- NULL
  }



  return(group_by_info)
}

#' @importFrom SummarizedExperiment colData

.extract_facet_by <- function(sce, facet_by) {
  if (!is.null(facet_by)) {
    if (!facet_by %in% colnames(SummarizedExperiment::colData(sce))) {
      stop("group_by is not a column name of sce's colData")
    }

    facet_by_info <- as.factor(colData(sce)[, facet_by])
    names(facet_by_info) <- colnames(sce)
  } else {
    facet_by_info <- NULL
  }
  return(facet_by_info)
}

.extract_filter_features <- function(exprsMat, feature_subset, n) {

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
  return(exprsMat)
}
