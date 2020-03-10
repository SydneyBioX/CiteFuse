
#' A function to perform DE analysis on CITE seq data
#'
#' @param sce A singlecellexperiment object
#' @param altExp_name A character indicates which expression matrix is used. by default is none (i.e. RNA).
#' @param exprs_value A character indicates which expression value in assayNames is used.
#' @param group A vector indicates the grouping of the data
#' @param method A character indicates the method used in DE analysis
#' @param exprs_pct A numeric indicates the threshold expression percentage of a gene to be considered in DE analysis
#' @param exprs_threshold A numeric indicates the threshold of expression. By default is 0.
#' @param return_all Whether return full list of DE genes
#' @param pval_adj A numeric indicates the threshold of adjusted p-value.
#' @param mean_diff A numeric indicates the threshold of difference of average expression.
#' @param pct_diff A numeric indicates the threshold of difference of percentage expression.
#' @param topN A numeric indicates the top number of genes will be included in the list.
#'
#' @importFrom randomForest randomForest
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom S4Vectors metadata
#'
#' @export
#'

DEgenes <- function(sce,
                    altExp_name = "none",
                    exprs_value = "logcounts",
                    group = NULL,
                    method = "wilcox",
                    exprs_pct = 0.1,
                    exprs_threshold = 0,
                    return_all = FALSE,
                    pval_adj = 0.05,
                    mean_diff = 0,
                    pct_diff = 0.1,
                    topN = 10) {

  method <- match.arg(method, c("wilcox"))

  if (is.null(group)) {
    stop("group is NULL.")
  }

  if (altExp_name != "none") {
    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
      stop("sce does not contain altExp_name as altExpNames")
    }

    if (!exprs_value %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
      stop("sce does not contain exprs_value as assayNames for altExp")
    }

    exprsMat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), exprs_value)

  } else {

    # if altExp_name is "none", then the assay in SingleCellExperiment is extracted (RNA in most of the cases)

    exprsMat <- SummarizedExperiment::assay(sce, exprs_value)
  }

  if (method == "wilcox") {
    de_res <- doWilcox(exprsMat,
                       cellTypes = group,
                       exprs_pct = exprs_pct,
                       exprs_threshold = exprs_threshold)
  }

  if (return_all) {
    return(de_res)
  } else {
    de_res <- selectDEgenes(de_res,
                            pval_adj = pval_adj,
                            mean_diff = mean_diff,
                            pct_diff = pct_diff,
                            topN = topN
    )

    return(de_res)
  }
}



#' A function to select DE genes
#'
#' @param de_res The de results returned by `DEgenes()``
#' @param pval_adj A numeric indicates the threshold of adjusted p-value.
#' @param mean_diff A numeric indicates the threshold of difference of average expression.
#' @param pct_diff A numeric indicates the threshold of difference of percentage expression.
#' @param topN A numeric indicates the top number of genes will be included in the list.
#'
#'
#' @export

selectDEgenes <- function(de_res,
                          pval_adj = 0.05,
                          mean_diff = 0,
                          pct_diff = 0.1,
                          topN = 10) {

  de_res <- lapply(de_res, function(x) {
    subset <- x[x$meanDiff > mean_diff &
                  x$p.adjust < pval_adj & x$pctDiff > pct_diff, ]
    subset <- subset[seq_len(min(nrow(subset), topN)),]
    subset
  })

  return(de_res)

}




doWilcox <- function(exprsMat, cellTypes, exprs_pct = 0.05, exprs_threshold = 0){
  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()

  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))


    meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
      rowMeans(exprsMat[, tmp_celltype == i])
    }))

    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      rowSums(exprsMat[, tmp_celltype == i] > exprs_threshold)/sum(tmp_celltype == i)
    }))

    meandiff <- meanPct[, 2] - meanPct[, 1]

    keep <- meanPct[,2] > exprs_pct

    test_res <- apply(exprsMat[keep, ], 1, function(x)
      stats::wilcox.test(x ~ tmp_celltype))

    test_res <- lapply(test_res, function(x) {
      stats = x$statistic
      pval = x$p.value
      p.adjust = p.adjust(pval, method = "BH")
      c(stats = stats, pval = pval, p.adjust = p.adjust)
    })

    # tt[[i]] <- topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    tt[[i]] <- do.call(rbind, test_res)

    rownames(tt[[i]]) <- rownames(exprsMat[keep, ])



    tt[[i]] <- cbind(tt[[i]], meanExprs.1 = meanExprs[rownames(tt[[i]]), 1])
    tt[[i]] <- cbind(tt[[i]], meanExprs.2 = meanExprs[rownames(tt[[i]]), 2])
    tt[[i]] <- cbind(tt[[i]], meanPct.1 = meanPct[rownames(tt[[i]]), 1])
    tt[[i]] <- cbind(tt[[i]], meanPct.2 = meanPct[rownames(tt[[i]]), 2])
    tt[[i]] <- cbind(tt[[i]], meanDiff = meandiff[rownames(tt[[i]])])
    tt[[i]] <- cbind(tt[[i]], pctDiff = meanPct[rownames(tt[[i]]), 2] -
                       meanPct[rownames(tt[[i]]), 1])

    tt[[i]] <- data.frame(tt[[i]])
    tt[[i]] <- tt[[i]][order(ifelse(tt[[i]]$meanDiff > 0, 0, 1), tt[[i]]$p.adjust),]
    tt[[i]] <- cbind(tt[[i]], name = rownames(tt[[i]]))
    tt[[i]] <- cbind(tt[[i]], group = levels(cellTypes)[i])
  }

  names(tt) <- levels(cellTypes)


  return(tt)
}



#' A function to generate circlepack plot to visualise the marker for each cluster
#'
#' @param de_list A list of results from `DE genes ()`
#'
#' @importFrom ggraph ggraph geom_node_circle geom_node_text geom_node_label
#' @import ggplot2
#'
#' @export

DEbubblePlot <- function(de_list) {

  type <- names(de_list)

  df_de <- do.call(rbind, lapply(seq_len(length(de_list)), function(i) {
    de_list[[i]] <- do.call(rbind, de_list[[i]])
    de_list[[i]]$type <- names(de_list)[i]
    de_list[[i]]
  }))

  p_value <- df_de$p.adjust
  p_value[p_value == 0] <- min(p_value[p_value > 0])

  data <- data.frame(
    root = rep("root", nrow(df_de)),
    group = paste("group", df_de$group, sep = "_"),
    subgroup = paste(df_de$group, df_de$type, sep = "_"),
    subsubgroup = df_de$name,
    type = df_de$type,
    value = rep(1, nrow(df_de)),
    size = -log10(p_value)
  )
  rownames(data) <- paste(data$subgroup, data$subsubgroup, sep = "|")


  edges <-  rbind(data.frame(from = data$root, to = data$group),
                  data.frame(from = data$group, to = paste(data$subgroup, data$subsubgroup, sep = "|")))

  labels <- unlist(lapply(strsplit(as.character(unique(unlist(edges))),"\\|"),
                          function(x) x[length(x)]))

  group <- unlist(lapply(strsplit(as.character(unique(unlist(edges))),"\\|"),
                         function(x) x[1]))

  vertices = data.frame(name = unique(unlist(edges)),
                        labels = labels,
                        group = group)

  fillColor <-  rep(NA, nrow(vertices))

  for (i in 1:length(type)) {
    fillColor[grep(type[i], vertices$group)] <- type[i]
  }

  fillColor[vertices$group == "root"] <- "0"
  fillColor[grep("group", vertices$group)] <- "group"

  vertices$fillColor <- fillColor

  vertices$size <- rep(1, nrow(vertices))
  rownames(vertices) <- vertices$name
  vertices[rownames(data), ]$size <- data$size


  vertices$originalGroup <- unlist(lapply(strsplit(as.character(vertices$group), "_"), "[[", 1))

  common <- rep(FALSE, nrow(vertices))
  group_list <- unique(vertices$originalGroup)
  group_list <- group_list[!group_list %in% c("root", "group")]

  for (i in 1:length(group_list)) {
    tmp <- vertices[vertices$originalGroup == group_list[i], ]
    common_features <- Reduce(intersect, lapply(lapply(unique(tmp$group), function(x) tmp[tmp$group == x, ]), function(x) x$labels))
    if (length(common_features) != 0) {
      common[vertices$originalGroup == group_list[i]][tmp$labels %in% common_features] <- TRUE
    }

  }


  vertices$common <- vertices$fillColor
  vertices$common <- ifelse(common, "common", vertices$common)
  vertices$fillColor <- factor(vertices$fillColor, levels = c("root", "group", type, "common"))

  mygraph <- igraph::graph_from_data_frame(edges, vertices = vertices )



  fill_colors <- c("white", "grey90", "#91B3D7", "#EA8783", "#ddb5d5")
  names(fill_colors) <- c("root", "group", type, "common")
  text_colors <- c("black", "black", "black", "black", "#D62728")
  names(text_colors) <- c("root", "group", type, "common")

  text_size <- vertices[!grepl("root|group", vertices$group), ]$size
  names(text_size) <- rownames(vertices[!grepl("root|group", vertices$group), ])
  text_size <- text_size/max(text_size) * 3
  text_size <- ifelse(text_size < 1.8, 1.8, text_size)

  g <- ggraph::ggraph(mygraph, layout = 'circlepack', weight = vertices$size) +
    ggraph::geom_node_circle(aes(fill = vertices$fillColor, color = vertices$common), size = 0.5) +
    # scale_fill_viridis_d() +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = text_colors) +
    theme_void() +
    ggraph::geom_node_text(aes(label = labels,
                               filter = leaf),
                           size = text_size) +
    ggraph::geom_node_label(aes(label = ifelse(grepl("group", as.character(vertices$labels)),
                                               as.character(vertices$labels), "")),
                            fill = "white",
                            repel = TRUE,
                            size = 3,
                            fontface = "bold") +
    theme(legend.position = "FALSE", aspect.ratio = 1)


  return(g)
}


#' A function to visualise the pairwise comparison of pvalue in different data modality.
#'
#' @param de_list A list including two lists results from `DE genes ()`.
#' @param feature_list A list including two lists features indicating the selected subset of features will be visualised
#'
#' @import ggplot2
#'
#' @export


DEcomparisonPlot <- function(de_list, feature_list) {

  if (length(de_list) != 2 | length(feature_list)  != 2) {
    stop("The length of de_list and feature_list not equal to 2")
  }

  #
  # for (i in 1:length(de_list)) {
  #   if (sum(!feature_list[[i]] %in% unique(as.character(unlist(lapply(de_list[[i]], "[[", "name"))))) != 0) {
  #     stop("Exist feature dost not have DE results.")
  #   }
  # }
  #
  de_pvalue_list <- list()

  for (i in 1:length(de_list)) {
    de_pvalue_list[[i]]  <- lapply(de_list[[i]], function(x) {
      p <- x[feature_list[[i]], ]$p.adjust
      p[is.na(p)] <- 1
      names(p) <- feature_list[[i]]
      p
    }
    )
  }

  df <- lapply(1:length(de_pvalue_list[[1]]), function(i) {

    df <- data.frame(log10(de_pvalue_list[[1]][[i]]),
                     -log10(de_pvalue_list[[2]][[i]]))
    df <- cbind(df, do.call(cbind, feature_list))
    colnames(df) <- c(names(de_list), paste(names(de_list), "name", sep = "_"))


    df <- df[rowSums(df[, 1:2]) != 0, ]

    df <- df[order(abs(df[, 2] - df[, 1]), decreasing = TRUE), ]

    df$group <- paste("Group", names(de_pvalue_list[[1]])[i])
    df
  })

  df <- do.call(rbind, df)

  g <- ggplot(data = df, aes(group = df$group)) +
    geom_segment(data = df, aes(x = 0, xend = df[, 1], y = 1:length(df[, 3]), yend = 1:length(df[, 3]))) +
    geom_segment(data = df, aes(x = 0, xend = df[, 2], y = 1:length(df[, 3]), yend = 1:length(df[, 3]))) +
    geom_point(data = df, aes(x = 0, y = 1:length(df[, 3])), color = "black", shape = 124, size = 10) +
    geom_point(data = df, aes(x = df[, 1], y = 1:length(df[, 3])), color = "red", size = 2) +
    geom_point(data = df, aes(x = df[, 2], y = 1:length(df[, 3])), color = "blue", size = 2) +
    scale_y_continuous(breaks = 1:(length(df[, 3])),
                       labels = c(as.character(df[, 3])),
                       sec.axis = sec_axis(~.,
                                           breaks = 1:(length(df[, 3])),
                                           labels = c(as.character(df[, 4])))) +
    facet_grid(group~., scales = "free_y", space = "free_y", switch = "both") +
    ylab("") + xlab("log10(p.adjust)                              -log10(p.adjust)") +  theme_bw() +
    theme(strip.placement = "outside", strip.text = element_text(size = 5)) +
    NULL
  return(g)

}





