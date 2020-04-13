#' ligandReceptorTest
#'
#' A function to perform ligand receptor analysis
#'
#'
#' @param sce A singlecellexperiment object
#' @param ligandReceptor_list A data.frame indicates the ligand receptor list
#' @param cluster A vector indicates the cluster results
#' @param RNA_exprs_value A character indicates which expression value
#' for RNA in assayNames is used.
#' @param use_alt_exp A logical vector indicates whether receptors
#' expression will use alternative expression matrix to quantify.
#' @param altExp_name A character indicates which expression matrix is used.
#' by default is ADT .
#' @param altExp_exprs_value A character indicates which expression value
#' in assayNames is used.
#' @param num_permute Number of permutation.
#' @param p_sig A numeric indicates threshold of the pvalue significance
#'
#' @return A SingleCellExperiment object with ligand receptor results
#'
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom S4Vectors metadata
#'
#' @examples
#' data(lr_pair_subset, package = "CiteFuse")
#' data(sce_control_subset, package = "CiteFuse")
#'
#' sce_control_subset <- normaliseExprs(sce = sce_control_subset,
#' altExp_name = "ADT",
#' transform = "zi_minMax")
#'
#' sce_control_subset <- normaliseExprs(sce = sce_control_subset,
#'                               altExp_name = "none",
#'                               exprs_value = "logcounts",
#'                               transform = "minMax")
#'
#' sce_control_subset <- ligandReceptorTest(sce = sce_control_subset,
#'                                   ligandReceptor_list = lr_pair_subset,
#'                                   cluster = sce_control_subset$SNF_W_louvain,
#'                                   RNA_exprs_value = "minMax",
#'                                   use_alt_exp = TRUE,
#'                                   altExp_name = "ADT",
#'                                   altExp_exprs_value = "zi_minMax",
#'                                   num_permute = 100)
#' @export


ligandReceptorTest <- function(sce,
                               ligandReceptor_list,
                               cluster,
                               RNA_exprs_value = "minMax",
                               use_alt_exp = TRUE,
                               altExp_name = "ADT",
                               altExp_exprs_value = "zi_minMax",
                               num_permute = 1000,
                               p_sig = 0.05) {


    if (!RNA_exprs_value %in% SummarizedExperiment::assayNames(sce)) {
        stop("sce does not contain RNA_exprs_value")
    }

    exprsMat1 <- SummarizedExperiment::assay(sce, RNA_exprs_value)

    if (use_alt_exp) {
        if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
            stop("sce does not contain altExp_name as altExpNames")
        }

        if (!altExp_exprs_value %in% assayNames(altExp(sce, altExp_name))) {
            stop("sce does not contain altExp_exprs_value as
                 assayNames for altExp")
        }

        # ADT exprssion matrix
        exprsMat2 <- assay(altExp(sce, altExp_name), altExp_exprs_value)

    } else {
        exprsMat2 <- assay(sce, RNA_exprs_value)
    }



    if (length(cluster) != ncol(sce)) {
        stop("The length of the cluster is not matched with
         the number of column of sce")
    }

    cluster_level <- levels(as.factor(droplevels(cluster)))

    keep <- ligandReceptor_list[, 1] %in% rownames(exprsMat1) &
        ligandReceptor_list[, 2] %in% rownames(exprsMat2)
    if (sum(keep) == 0) {
        stop("None of the ligand-receptor pairs are in the provided sce")
    }
    ligandReceptor_list <- ligandReceptor_list[keep, ]


    l_zeros <- vapply(cluster_level, function(x)
        apply(exprsMat1[ligandReceptor_list[, 1], cluster == x], 1, .zeroProp),
        numeric(nrow(exprsMat1[ligandReceptor_list[, 1], ])))
    r_zeros <- vapply(cluster_level, function(x)
        apply(exprsMat2[ligandReceptor_list[, 2], cluster == x], 1, .zeroProp),
        numeric(nrow(exprsMat2[ligandReceptor_list[, 2], ])))



    l_express <- l_zeros < 0.9
    r_express <- r_zeros < 0.9


    keep_lr <- rowSums(l_express) != 0 & rowSums(r_express) != 0

    ligandReceptor_list <- ligandReceptor_list[keep_lr, ]


    lr_mean <- lapply(seq_len(nrow(ligandReceptor_list)), function(p) {


        pair <- ligandReceptor_list[p, ]
        res <- vapply(cluster_level, function(x) {
            rna <- mean((exprsMat1[pair[1], as.character(cluster) == x]))
            adt <- mean((exprsMat2[pair[2], as.character(cluster) == x]))
            df <- c(rna, adt)

        }, numeric(2))
        rownames(res) <- pair
        res
    })



    permute_lr_mean <- list()




    permute_lr_mean <- lapply(seq_len(num_permute), function(idx_per) {

        if (idx_per %% 100 == 0) cat(idx_per, "......")

        cluster_permute <- cluster[sample(length(cluster), length(cluster))]
        l <- lapply(seq_len(nrow(ligandReceptor_list)), function(p) {

            pair <- ligandReceptor_list[p, ]
            res <- vapply(cluster_level, function(x) {
                rna <- mean((exprsMat1[pair[1],
                                       as.character(cluster_permute) == x]))
                adt <- mean((exprsMat2[pair[2],
                                       as.character(cluster_permute) == x]))
                df <- c(rna, adt)

            }, numeric(2))
            rownames(res) <- pair
            res
        })
        l
    })




    pvalue <- list()

    for (pair in seq_len(nrow(ligandReceptor_list))) {

        res <- lapply(permute_lr_mean, function(x) x[[pair]])
        res <- lapply(res, function(x) .get_expand_grid_average(x[1, ], x[2, ]))
        res <- do.call(rbind, res)

        observed <- .get_expand_grid_average(lr_mean[[pair]][1, ],
                                             lr_mean[[pair]][2, ])

        pvalue[[pair]] <- vapply(seq_len(length(observed)), function(idx) {
            p <- 1 - sum(observed[idx] > res[, idx])/length(res[, idx])
        }, numeric(1))
    }



    pvalue <- do.call(rbind, pvalue)

    rownames(pvalue) <- paste(ligandReceptor_list[, 1],
                              ligandReceptor_list[, 2], sep = "|")

    pvalue_filter <- pvalue[apply(pvalue, 1, function(x) sum(x < p_sig)) != 0,]


    grid_list <- expand.grid(seq_len(length(cluster_level)),
                             seq_len(length(cluster_level)))

    colnames(pvalue_filter) <- paste(grid_list[, 1],
                                     grid_list[, 2], sep = "|")

    pvalue_filter[pvalue_filter == 0] <- 1/(10*num_permute)




    num_interaction <- colSums(pvalue_filter < p_sig)

    df_cci <- do.call(rbind, strsplit(names(num_interaction), "\\|"))
    df_cci <- cbind(df_cci, num_interaction)

    df_cci <- apply(df_cci, 2, as.numeric)

    mat_cci <- matrix(0, length(cluster_level),  length(cluster_level))
    for (i in seq_len(nrow(mat_cci))) {
        for (j in seq_len(ncol(mat_cci))) {
            mat_cci[i, j] <- df_cci[df_cci[, 1] == i & df_cci[, 2] == j, 3]
        }
    }

    colnames(mat_cci) <- rownames(mat_cci) <-
        paste("Group", seq_len(length(cluster_level)), sep = "_")

    mat_cci_sym <- mat_cci + t(mat_cci)
    diag(mat_cci_sym) <- diag(mat_cci)


    idx_sig <- which(pvalue_filter < p_sig, arr.ind = TRUE)
    df_sig_list <- data.frame(lr_pair = rownames(idx_sig),
                              cluster_pair =
                                  colnames(pvalue_filter)[idx_sig[, 2]])

    df_sig_list$ligand <- unlist(lapply(strsplit(rownames(idx_sig),
                                                 "\\|"), "[[", 1))
    df_sig_list$receptor <- unlist(lapply(strsplit(rownames(idx_sig),
                                                   "\\|"), "[[", 2))
    ligand_cluster <- unlist(lapply(strsplit(
        as.character(df_sig_list$cluster_pair),
        "\\|"), "[[", 1))

    df_sig_list$ligand_cluster <- ligand_cluster

    receptor_cluster <- unlist(lapply(strsplit(
        as.character(df_sig_list$cluster_pair),
        "\\|"), "[[", 2))

    df_sig_list$receptor_cluster <- receptor_cluster



    LRanalysis_results <- list()
    LRanalysis_results$LRanalysis_pvalue <- pvalue_filter
    LRanalysis_results$LRanalysis_group <- cluster
    LRanalysis_results$LRanalysis_pairsCount_sym <- mat_cci_sym
    LRanalysis_results$LRanalysis_pairsCount <- mat_cci
    LRanalysis_results$LRanalysis_pairsList <- df_sig_list
    LRanalysis_results$receptor_type <- ifelse(use_alt_exp,
                                               altExp_name,
                                               "RNA")

    res_meta_name <- paste("LRanalysis_results",
                           LRanalysis_results$receptor_type,
                           sep = "_")

    S4Vectors::metadata(sce)[[res_meta_name]] <- LRanalysis_results

    return(sce)

}



.zeroProp <- function(x) {
    sum(x == 0)/length(x)
}







.get_expand_grid_average <- function(x, y) {
    grid_list <- expand.grid(seq_len(length(x)), seq_len(length(x)))
    res <- (x[grid_list[, 1]] + y[grid_list[, 2]])/2
    names(res) <- paste(grid_list[, 1], grid_list[, 2], sep = "|")


    return(res)
}

#' visLigandReceptor
#'
#' A function to visualise ligand receptor analysis
#'
#'
#' @param sce A singlecellexperiment object
#' @param type A character indicates the type of the plot for ligand receptor
#' restuls visualisation, option includes "pval_heatmap", "pval_dotplot",
#' "group_network", "group_heatmap", and "lr_network"
#' @param receptor_type A character indicates which receptor expression's
#' ligand receptor results are used to generate the figures.
#'
#' @return A plot visualise the ligand receptor results
#'
#'
#' @importFrom S4Vectors metadata
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 melt
#' @importFrom pheatmap pheatmap
#' @importFrom igraph graph_from_data_frame V E layout_in_circle
#' @importFrom graphics plot legend
#' @import ggplot2
#'
#' @examples
#' data(lr_pair_subset, package = "CiteFuse")
#' data(sce_control_subset, package = "CiteFuse")
#'
#' sce_control_subset <- normaliseExprs(sce = sce_control_subset,
#' altExp_name = "ADT",
#' transform = "zi_minMax")
#'
#' sce_control_subset <- normaliseExprs(sce = sce_control_subset,
#'                               altExp_name = "none",
#'                               exprs_value = "logcounts",
#'                               transform = "minMax")
#'
#' sce_control_subset <- ligandReceptorTest(sce = sce_control_subset,
#'                                   ligandReceptor_list = lr_pair_subset,
#'                                   cluster = sce_control_subset$SNF_W_louvain,
#'                                   RNA_exprs_value = "minMax",
#'                                   use_alt_exp = TRUE,
#'                                   altExp_name = "ADT",
#'                                   altExp_exprs_value = "zi_minMax",
#'                                   num_permute = 100)
#' visLigandReceptor(sce_control_subset,
#' type = "pval_heatmap",
#' receptor_type = "ADT")
#' @export

visLigandReceptor <- function(sce,
                              type = c("pval_heatmap", "pval_dotplot",
                                       "group_network", "group_heatmap",
                                       "lr_network"),
                              receptor_type = NULL) {


    type <- match.arg(type, c("pval_heatmap", "pval_dotplot", "group_network",
                              "group_heatmap", "lr_network"))

    if (is.null(receptor_type)) {
        receptor_type <- "ADT"
    }

    meta_name <- paste("LRanalysis_results",
                       receptor_type,
                       sep = "_")

    if (!meta_name %in% names(S4Vectors::metadata(sce))) {
        stop("Please perform ligandReceptorTest() first")
    }


    LRanalysis_results <- S4Vectors::metadata(sce)[[meta_name]]

    pvalue_filter <- as.matrix(LRanalysis_results$LRanalysis_pvalue)

    cluster_level <- levels(factor(droplevels(
        LRanalysis_results$LRanalysis_group)))



    cluster_colors <- cite_colorPal(length(cluster_level))
    names(cluster_colors) <- cluster_level

    if (type == "pval_heatmap") {

        reds <- c("#67000D", "#A50F15", "#CB181D",
                  "#EF3B2C", "#FB6A4A", "#FC9272",
                  "#FCBBA1", "#FEE0D2", "#FFF5F0")

        pval_colors <- c(grDevices::colorRampPalette(reds[seq_len(4)])(100),
                         grDevices::colorRampPalette(reds[seq(4, 9)])(400),
                         rep("#FFF5F0", 500))

        ligand_cluster <-
            factor(unlist(lapply(strsplit(colnames(pvalue_filter),
                                          "\\|"), "[[", 1)),
                   levels = cluster_level)

        receptor_cluster <-
            factor(unlist(lapply(strsplit(colnames(pvalue_filter),
                                          "\\|"), "[[", 2)),
                   levels = cluster_level)

        annotation_row <- data.frame(ligand_cluster = ligand_cluster,
                                     receptor_cluster = receptor_cluster)



        rownames(annotation_row) <- colnames(pvalue_filter)

        annotation_colors <- list(ligand_cluster = cluster_colors,
                                  receptor_cluster = cluster_colors)



        pheatmap::pheatmap(t(pvalue_filter),
                           color = pval_colors,
                           annotation_row = annotation_row,
                           annotation_colors = annotation_colors,
                           border_color = "grey100")

    }

    if (type == "pval_dotplot") {

        df_toPlot <- reshape2::melt(pvalue_filter)

        g <- ggplot(df_toPlot, aes(x = df_toPlot$Var1, y = df_toPlot$Var2)) +
            geom_point(aes(size = -log10(df_toPlot$value),
                           color = -log10(df_toPlot$value))) +
            scale_color_viridis_c(direction = 1, end = 0.95) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1.1),
                  axis.text.y = element_text(size = 8),
                  axis.ticks = element_blank(),
                  legend.position = "bottom") +
            labs(size = "-log10(p)", color = "-log10(p)") +
            ylab("Cluster Pairs") + xlab("Ligand Receptor Pairs")
        return(g)
    }


    if (type == "lr_network") {
        df_sig_list <- LRanalysis_results$LRanalysis_pairsList

        df_sig_list[, 3] <- paste("L", df_sig_list[, 3], sep = "_")
        df_sig_list[, 4] <- paste("R", df_sig_list[, 4], sep = "_")
        g <- igraph::graph_from_data_frame(df_sig_list[, 3:4],
                                           directed = FALSE)



        igraph::V(g)$label <- unlist(lapply(strsplit(names(V(g)), "_"),
                                            function(x) paste(x[-1],
                                                              collapse = "_")))
        igraph::V(g)$class <- unlist(lapply(strsplit(names(V(g)), "_"),
                                            "[[", 1))
        numeric_class <- as.numeric(as.factor(V(g)$class))

        igraph::V(g)$type <- c(TRUE, FALSE)[numeric_class]
        igraph::V(g)$shape <- c("circle", "square")[numeric_class]
        igraph::V(g)$color <- c("#A0CBE8", "#FFBE7D")[numeric_class]
        igraph::V(g)$size <- 5
        igraph::V(g)$label.cex <- 0.4
        igraph::V(g)$label.color <- "black"
        df_sig_list$cluster_pair <- as.factor(df_sig_list$cluster_pair)
        nlevels_cluster_pair <- nlevels(df_sig_list$cluster_pair)
        num_cluster_pair <- as.numeric(df_sig_list$cluster_pair)
        igraph::E(g)$color <-
            cite_colorPal(nlevels_cluster_pair)[num_cluster_pair]


        graphics::plot(g)

        graphics::legend('topleft',
                         legend = c(levels(df_sig_list$cluster_pair)),
                         col = cite_colorPal(nlevels_cluster_pair),
                         pch = 95,
                         ncol = ceiling(nlevels_cluster_pair/10),
                         cex = 0.5)


    }

    if (type == "group_network") {



        mat_cci <- LRanalysis_results$LRanalysis_pairsCount_sym

        mat_cci <- reshape2::melt(mat_cci)

        mat_cci <- mat_cci[mat_cci$value != 0, ]


        g <- igraph::graph_from_data_frame(mat_cci)

        tab <- as.numeric(mat_cci$value)
        l <- igraph::layout_in_circle(g)
        igraph::E(g)$width <- (tab/sum(tab))

        igraph::E(g)$width <- igraph::E(g)$width/max(igraph::E(g)$width) * 10
        igraph::E(g)$arrow.size <- 0.5
        igraph::V(g)$color <- cluster_colors[gsub("Group_", "", names(V(g)))]

        plot(g, layout = l)
    }

    if (type == "group_heatmap") {

        mat_cci <- LRanalysis_results$LRanalysis_pairsCount
        pheatmap::pheatmap(mat_cci,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           border_color = "grey100")
    }

}
