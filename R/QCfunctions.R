#' A function to preprocess the list of expression matrix
#'
#' @description  This function will keep the samples that are common across the list of expression matrix,
#' and filter the features that are all zeros across samples, and finally construct a \code{SingleCellExperiment}
#' object
#'
#' @param exprsMat A list or a matrix indicates the expression matrices of the
#' testing datasets (each matrix must be \code{matrix} or \code{dgCMatrix} class)
#' @param return_sce A logical input indicates whether a \code{SingleCellExperiment}
#' object will be return
#' @param assay_matrix A integer indicates which list will be used as `assay` input of `SingleCellExperiment`
#' @param filter_features A logical input indicates whether the features with all zeros will be removed
#' @param rowData A DataFrame indicates the rowData to be stored in the sce object
#' @param colData A DataFrame indicates the colData to be stored in the sce object
#'
#' @return either a SingleCellExperiment object or a preprocessed expression matrix
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix rowSums
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods as is
#'
#' @export



preprocessing <- function(exprsMat = NULL,
                          return_sce = TRUE,
                          assay_matrix = 1,
                          filter_features = TRUE,
                          rowData = NULL,
                          colData = NULL) {


  if (!any("matrix" %in% methods::is(exprsMat),
           "dgCMatrix" %in% methods::is(exprsMat),
           "list" %in% methods::is(exprsMat))) {
    stop("exprsMat need to be a matrix or a list")
  }

  if ("list" %in% methods::is(exprsMat)) {
    if (any(!unlist(lapply(exprsMat, function(x) "dgCMatrix" %in% is(x))) &
            !unlist(lapply(exprsMat, function(x) "matrix" %in% is(x))))) {
      stop("Please make sure every expression matrix in
           the list are matrix or sparse matrix")
    }
  }

  if (any("matrix" %in% methods::is(exprsMat),
          "dgCMatrix" %in% methods::is(exprsMat))) {
    exprsMat <- list(exprsMat)
  }

  common_cells <- Reduce(intersect, lapply(exprsMat, colnames))

  if (length(common_cells) == 0) {
    stop("There is no common cells in this list... please check the matrix input")
  }




  # only keep the samples that are common across the list of exprsMat
  exprsMat <- lapply(exprsMat, function(exprs) {
    exprs <- exprs[, common_cells]
    if (!"dgCMatrix" %in% methods::is(exprs)) {
      exprs <- methods::as(exprs, "dgCMatrix")
    }
    exprs
  })

  # filter the features with all zeros

  if (filter_features) {
    exprsMat <- lapply(exprsMat, function(exprs) {
      if ("dgCMatrix" %in% methods::is(exprsMat)) {
        rowsums <- Matrix::rowSums(exprs)
        exprs <- exprs[rowsums != 0, ]

      }
      exprs
    })
  }



  if (return_sce) {

    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = exprsMat[[assay_matrix]]))

    list_idx <- seq_len(length(exprsMat))[-assay_matrix]

    for (i in list_idx) {
      if (is.null(names(exprsMat)[i])) {
        name_exprs <- paste("altExp", i, sep = "_")
      } else {
        name_exprs <- names(exprsMat)[i]
      }
      SingleCellExperiment::altExp(sce, name_exprs) <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = exprsMat[[i]]))
    }

    if (!is.null(rowData)) {

      if (!all(rownames(sce) %in% rownames(rowData))) {
        stop("Some rownames of the assay matrix does not have info in rowData (rownames of rowData)")
      }

      SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(rowData[rownames(sce), ])
    }

    if (!is.null(colData)) {

      if (!all(colnames(sce) %in% rownames(colData))) {
        stop("Some colnames of the assay matrix does not have info in colData (rownames of colData)")
      }

      SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(colData[colnames(sce), ])
    }

    return(sce)

  } else {
    return(exprsMat)
  }

}



#' readFrom10X
#'
#' A function to read the data from 10X
#'
#'
#' @param dir A character indicates the directory of the 10X files
#' @param type A character indicates the format of the data, sparse or HDF5
#' @param feature_named_by A character indicates whehter the genes will be named by gene_id or gene_symbol
#' @param filter_features A logical input indicates whether the features with all zeros will be removed
#'
#' @return a SingleCellExperiment object
#'
#' @importFrom rhdf5 h5read
#' @importFrom Matrix readMM sparseMatrix
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods as
#' @importFrom utils read.delim
#'
#' @export
#'


readFrom10X <- function(dir,
                        type = c("auto", "sparse", "HDF5"),
                        feature_named_by = c("gene_id", "gene_symbol"),
                        filter_features = TRUE) {





  type <- match.arg(type, c("auto", "sparse", "HDF5"))
  feature_named_by <- match.arg(feature_named_by, c("gene_id", "gene_symbol"))

  all_files <- list.files(dir)


  if (type == "auto") {

    if (sum(grepl(c("barcodes.tsv|features.tsv|matrix.mtx"), all_files)) == 3) {
      type <- "sparse"
    } else if (any(grepl("\\.h5", all_files))) {
      type <- "HDF5"
    } else {
      stop("No avaliable files are found in the given dir.")
    }

  }

  if (type == "sparse") {
    if (sum(grepl(c("barcodes.tsv|features.tsv|matrix.mtx"), all_files)) != 3) {
      stop("No avaliable files are found in the given dir.")
    }
  }


  if (type == "HDF5") {
    if (!any(grepl("\\.h5", all_files))) {
      stop("No avaliable files are found in the given dir.")
    }
  }



  if (type == "HDF5") {

    h5_path <- paste0(dir, list.files(dir, pattern = "\\.h5"))

    #rhdf5::h5ls(h5_path, all = TRUE)

    barcodes  <- as.character(h5read(h5_path, paste0("matrix", "/barcodes")))


    gene_id  <- as.character(h5read(h5_path, paste0("matrix", "/features/id")))
    gene_symbol  <- as.character(h5read(h5_path, paste0("matrix", "/features/name")))
    gene_type  <- as.character(h5read(h5_path, paste0("matrix", "/features/feature_type")))

    features <- data.frame(gene_id = gene_id,
                           gene_symbol = gene_symbol,
                           gene_type = gene_type)

    counts <- rhdf5::h5read(h5_path, paste0("matrix", "/data"))
    indices <- rhdf5::h5read(h5_path, paste0("matrix", "/indices"))
    indptr <- rhdf5::h5read(h5_path, paste0("matrix", "/indptr"))
    shape <- rhdf5::h5read(h5_path, paste0("matrix", "/shape"))

    exprs_mat <- Matrix::sparseMatrix(i = indices + 1,
                                      p = indptr,
                                      x = as.numeric(x = counts),
                                      dims = shape,
                                      giveCsparse = TRUE)

  }


  if (type == "sparse") {

    mat_path <- paste0(dir, "matrix.mtx.gz")
    exprs_mat <- Matrix::readMM(mat_path)
    exprs_mat <- methods::as(exprs_mat, "dgCMatrix")

    barcodes_path <- paste0(dir, "barcodes.tsv.gz")
    barcodes <- utils::read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
    barcodes <- barcodes[, 1]



    feature_path <- paste0(dir, "features.tsv.gz")
    features <- utils::read.delim(feature_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(features) <- c("gene_id", "gene_symbol", "gene_type")




  }

  rownames(exprs_mat) <- features[, feature_named_by]
  colnames(exprs_mat) <- barcodes


  gene_idx <- features$gene_type == "Gene Expression"
  adt_idx <- features$gene_type == "Antibody Capture"

  exprs_rna <- exprs_mat[gene_idx, ]
  exprs_adt <- exprs_mat[adt_idx, ]

  rownames(features) <- features$gene_id


  sce <- preprocessing(list(RNA = exprs_rna, ADT = exprs_adt),
                       rowData = features[gene_idx, ],
                       filter_features = filter_features)




  return(sce)

}




#' normaliseExprs
#'
#' A function that perform normalisation for alternative expression
#'
#' @param sce A \code{SingleCellExperiment} object
#' @param altExp_name Name of alternative expression that will be used to perform normalisation
#' @param exprs_value A character indicates which expression value in assayNames is used.
#' @param transform type of transformation, either log or clr (Centered log ratio transform)
#' @param log_offset Numeric scalar specifying the pseudo-count to add when log-transforming expression values. Default is 1
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom Matrix rowMeans
#'
#' @return a SingleCellExperiment object
#'
#' @export

normaliseExprs <- function(sce,
                           altExp_name = NULL,
                           exprs_value = "counts",
                           transform = c("log", "clr", "zi_minMax", "minMax"),
                           log_offset = NULL

) {

  transform <- match.arg(transform, c("log", "clr", "zi_minMax", "minMax"),
                         several.ok = TRUE)

  if (altExp_name != "none") {
    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
      stop("sce does not contain altExp_name as altExpNames")
    }

    assaynames <- SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))
    if (!exprs_value %in% assaynames) {
      stop("sce does not contain exprs_value as assayNames for altExp")
    }

    exprs <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), exprs_value)

  } else {

    # if altExp_name is "none", then the assay in SingleCellExperiment is extracted (RNA in most of the cases)

    assaynames <- SummarizedExperiment::assayNames(sce)
    if (!exprs_value %in% assaynames) {
      stop("sce does not contain exprs_value as assayNames")
    }

    exprs <- SummarizedExperiment::assay(sce, exprs_value)
  }

  if (is.null(log_offset)) {
    log_offset <- 1
  }

  if (transform %in% "log") {

    exprs_norm <- log(exprs + log_offset)

  }

  if (transform %in% "clr") {

    exprs_norm <- .clr(exprs)

  }

  if (transform %in% "zi_minMax") {

    # if (!"logcounts" %in% SummarizedExperiment::assayNames(SingleCellExperiment::altExp(sce, altExp_name))) {
    #   exprs_log <- log(exprs + log_offset)
    # } else {
    #   exprs_log <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), "logcounts")
    # }
    #
    exprs_norm <- apply(exprs, 1, .ziMinMax)
    exprs_norm <- t(exprs_norm)

  }

  if (transform %in% "minMax") {

    # if (!"logcounts" %in% assaynames) {
    #   exprs_log <- log(exprs + log_offset)
    # } else {
    #   exprs_log <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), "logcounts")
    # }

    exprs_norm <- apply(exprs, 1, .minMax)
    exprs_norm <- t(exprs_norm)

  }


  if (transform == "log") {
    new_assay_name <- "logcounts"
  } else {
    new_assay_name <- transform
  }


  if (altExp_name != "none") {
    SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name),
                                new_assay_name) <- exprs_norm
  } else {
    SummarizedExperiment::assay(sce, new_assay_name) <- exprs_norm
  }

  return(sce)
}


.minMax <- function(x) {
  if (max(x) != 0) {
    res <- (x - min(x))/(max(x) - min(x))
  } else {
    res <- rep(0, length(x))
  }
  return(res)
}

.ziMinMax <- function(x) {
  res <- .minMax(x)
  res <- ifelse(res < 0.5, 0, res)
  return(res)
}

.clr <- function(X) {

  X[X == 0] <- min(X[X != 0])
  logX <- log(X)
  logSet <- logX[, seq_len(ncol(X)), drop = FALSE]
  ref <- Matrix::rowMeans(logSet)
  res <- sweep(logX, 1, ref, "-")

  return(res)
}


#' crossSampleDoublets
#'
#' A function that perform normalisation for alternative expression
#'
#' @param sce A \code{SingleCellExperiment} object
#' @param altExp_name Name of alternative expression that will be used to perform normalisation.
#' If it is NULL, it will set to HTO.
#' @param totalExp_threshold the threshold indicates for the HTO less than this threshold
#' will be filtered from the analysis
#'
#' @return A SingleCellExperiment Object
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom Matrix rowSums
#' @importFrom mixtools normalmixEM
#' @importFrom S4Vectors metadata
#' @importFrom stats kmeans
#' @importFrom methods is
#'
#'
#' @export


crossSampleDoublets <- function(sce,
                                altExp_name = NULL,
                                totalExp_threshold = 10) {

  if (is.null(altExp_name)) {
    if (!"HTO" %in% SingleCellExperiment::altExpNames(sce)) {
      stop("There is no HTO data in the object")
    } else {
      altExp_name <- "HTO"
    }
  }

  if (!"logcounts" %in% SummarizedExperiment::assayNames(altExp(sce, altExp_name))) {
    warning("HTO does not contain logcounts... we will perform normaliseExprs() to get logcounts")
    sce <- normaliseExprs(sce, altExp_name, "log")
  }

  hto_cellHash_log <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), "logcounts")

  hto_cellHash_log <- hto_cellHash_log[Matrix::rowSums(hto_cellHash_log) >
                                         totalExp_threshold, ]

  hto_threshold <- list()
  for (i in seq_len(nrow(hto_cellHash_log))) {

    vec <- hto_cellHash_log[i,][hto_cellHash_log[i,] > 0]
    mixmdl <- try(mixtools::normalmixEM(vec,
                                        fast = TRUE, maxrestarts = 1000,
                                        k = 2, maxit = 10000,
                                        mu = c(0, 10),
                                        lambda = c(1/2),
                                        sigma = rep(2, 2),
                                        ECM = TRUE, verb = FALSE),
                  silent = TRUE)



    if ("try-error" %in% methods::is(mixmdl)) {
      km <- stats::kmeans(vec, centers = 2)
      hto_threshold[[i]] <- min(max(vec[km$cluster == 1]), max(vec[km$cluster == 2]))
    } else {
      hto_threshold[[i]] <- getThreshold(mixmdl)
    }
  }

  names(hto_threshold) <- NULL
  hto_threshold <- unlist(hto_threshold)


  hto_cellHash_pass <- sapply(seq_len(nrow(hto_cellHash_log)), function(x) {
    hto_cellHash_log[x,] > hto_threshold[x]
  })


  colnames(hto_cellHash_pass) <- rownames(hto_cellHash_log)


  hto_cellHash_mix_label <- apply(hto_cellHash_pass, 1, function(x) {
    if (sum(x) == 0) {
      "negative"
    } else if (sum(x) == 1) {
      which(x)
    } else {
      "doublet/multiplet"
    }
  })

  doubletClassify_between_class <- ifelse(!hto_cellHash_mix_label %in%
                                            c("negative", "doublet/multiplet"),
                                          "Singlet", hto_cellHash_mix_label)


  sce$doubletClassify_between_label <- hto_cellHash_mix_label
  sce$doubletClassify_between_class <- doubletClassify_between_class

  S4Vectors::metadata(sce) <- c(S4Vectors::metadata(sce),
                                list(doubletClassify_between_threshold = hto_threshold,
                                     doubletClassify_between_resultsMat = hto_cellHash_pass))


  return(sce)

}



#' plotHTO
#'
#' A function to plot HTO expression
#'
#' @param sce sce
#' @param which_idx which_idx
#' @param altExp_name altExp_name
#' @param ncol ncol
#'
#' @return A plot visualising the HTO expression
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom gridExtra grid.arrange
#' @importFrom cowplot axis_canvas insert_xaxis_grob insert_yaxis_grob
#' @importFrom grid unit
#' @importFrom utils combn
#' @import ggplot2
#'
#' @export


plotHTO <- function(sce,
                    which_idx = seq_len(2),
                    altExp_name = NULL,
                    ncol = 2) {

  combination <- utils::combn(which_idx, 2)

  ggList <- apply(combination, 2, function(x) {
    plotHTOSingle(sce, which_idx = x, altExp_name = altExp_name)
  })



  do.call(gridExtra::grid.arrange, c(ggList, ncol = min(length(ggList), ncol)))
}


#' plotHTOSingle
#'
#' A function to plot HTO expression
#'
#' @param sce sce
#' @param which_idx which_idx
#' @param altExp_name altExp_name
#'
#' @return A plot visualising the HTO expression
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom cowplot axis_canvas insert_xaxis_grob insert_yaxis_grob
#' @importFrom grid unit




plotHTOSingle <- function(sce,
                          which_idx = seq_len(2),
                          altExp_name = NULL
) {

  if (is.null(altExp_name)) {
    if (!"HTO" %in% altExpNames(sce)) {
      stop("There is no HTO data in the object")
    } else {
      altExp_name <- "HTO"
    }
  }

  if (!"logcounts" %in% SummarizedExperiment::assayNames(altExp(sce, altExp_name))) {
    warning("HTO does not contain logcounts... we will perform normaliseExprs() to get logcounts")
    sce <- normaliseExprs(sce, altExp_name, "log")
  }

  noThreshold <- FALSE
  if (!"doubletClassify_between_label" %in% colnames(SummarizedExperiment::colData(sce))) {
    warning("Haven't performed doubletClassify() yet!")
    noThreshold <- TRUE
  }

  hto_cellHash_log <- assay(SingleCellExperiment::altExp(sce, altExp_name), "logcounts")


  df <- data.frame(t(as.matrix(hto_cellHash_log[which_idx, ])))

  colnames(df) <- lapply(strsplit(rownames(hto_cellHash_log[which_idx, ]), "\\."), "[[", 1)


  if (!noThreshold) {
    hto_cellHash_pass <- metadata(sce)[["doubletClassify_between_resultsMat"]]

    doublets <- hto_cellHash_pass[, which_idx[1]] & hto_cellHash_pass[, which_idx[2]]

    hto_threshold <- metadata(sce)[["doubletClassify_between_threshold"]]

  } else {
    doublets <- rep(FALSE, ncol(sce))
    hto_threshold <- rep(NULL, 2)
    hto_cellHash_pass <- NULL
  }


  pmain <- ggplot(df, aes(x = df[, 1], y = df[, 2], color = doublets)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = hto_threshold[which_idx[2]], col = "red",
               linetype = 2, size = 1) +
    geom_vline(xintercept = hto_threshold[which_idx[1]], col = "red",
               linetype = 2, size = 1) +
    scale_color_manual(values = c("#377EB8", "#E41A1C")) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2])




  xdens <- cowplot::axis_canvas(pmain, axis = "x") +
    geom_density(data = df, aes(x = df[, 1],
                                fill = hto_cellHash_pass[, which_idx[1]],
                                alpha = 0.5)) +
    scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
    NULL

  ydens <- cowplot::axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
    geom_density(data = df, aes(x = df[, 2],
                                fill = hto_cellHash_pass[, which_idx[2]],
                                alpha = 0.5)) +
    coord_flip() +
    scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
    NULL

  p1 <- cowplot::insert_xaxis_grob(pmain, xdens,
                                   grid::unit(.2, "null"), position = "top")

  p2 <- cowplot::insert_yaxis_grob(p1, ydens,
                                   grid::unit(.2, "null"), position = "right")

  # ggdraw(p2)
  return(p2)
}




#' withinSampleDoublets
#'
#' doublet identification within batch
#'
#' @param sce a SingleCellExperiment
#' @param altExp_name expression name of HTO matrix
#' @param eps eps of DBSCAN
#' @param minPts minPts of DBSCAN
#'
#' @return A SingleCellExperiment object
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment altExpNames altExp counts
#' @importFrom dbscan dbscan
#' @export


withinSampleDoublets <- function(sce,
                                 altExp_name = NULL,
                                 eps = 200,
                                 minPts = 50) {
  sce$nUMI <- Matrix::colSums(SingleCellExperiment::counts(sce))

  if (is.null(altExp_name)) {
    if (!"HTO" %in% SingleCellExperiment::altExpNames(sce)) {
      stop("There is no HTO data in the object")
    } else {
      altExp_name <- "HTO"
    }
  }


  if (!"logcounts" %in% SummarizedExperiment::assayNames(altExp(sce, altExp_name))) {
    warning("HTO does not contain logcounts...
            we will perform normaliseExprs() to get logcounts")
    sce <- normaliseExprs(sce, altExp_name, "log")
  }

  if (!"doubletClassify_between_label" %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("Haven't performed doubletClassify_between() yet!")
  }

  hto_threshold <- metadata(sce)[["doubletClassify_between_threshold"]]


  hto_cellHash_log <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), "logcounts")

  hto_cellHash_mix_label <- sce$doubletClassify_between_label

  nUMI <- sce$nUMI

  batch_doublets_list <- lapply(seq_len(nrow(hto_cellHash_log)), function(i) {

    db_cluster <- dbscan::dbscan(cbind(nUMI, hto_cellHash_log[i,]),
                                 eps = eps, minPts = minPts)
    batch_doublets <- (db_cluster$cluster == 0 &
                         hto_cellHash_log[i,] > hto_threshold[i] &
                         hto_cellHash_mix_label != "doublet/multiplet")

  })

  batch_doublets_mat <- t(do.call(rbind, batch_doublets_list))


  doubletClassify_within_label <- apply(batch_doublets_mat, 1, function(res) {
    if (sum(res, na.rm = TRUE) == 0) {
      "NotDoublets(Within)"
    } else {
      paste("Doublets(Within)", which(res), sep = "_")
    }
  })

  doubletClassify_within_class <- ifelse(doubletClassify_within_label == "NotDoublets(Within)", "Singlet", "Doublet")

  sce$doubletClassify_within_label <- doubletClassify_within_label
  sce$doubletClassify_within_class <- doubletClassify_within_class

  S4Vectors::metadata(sce) <- c(S4Vectors::metadata(sce),
                                list(doubletClassify_within_resultsMat = batch_doublets_mat))

  return(sce)

}


#' @importFrom stats uniroot
#' @importFrom methods is


getThreshold <- function(mixmdl, verbose = FALSE){

  # if (verbose) {
  #   plot(mixmdl, which = 2)
  # }

  membership <- apply(mixmdl$posterior, 1, which.max)
  m_list <- sort(unique(membership))

  mu_list <- mixmdl$mu
  names(mu_list) <- seq_len(length(mu_list))
  mu_list <- mu_list[m_list]

  if (length(mu_list) > 1) {
    idx1 <- as.numeric(names(mu_list)[order(mu_list)][1])
    idx2 <- as.numeric(names(mu_list)[order(mu_list)][2])

    root <- try(stats::uniroot(funMixModel,
                               interval = c(mixmdl$mu[idx1] - mixmdl$sigma[idx1], mixmdl$mu[idx2] + mixmdl$sigma[idx2]),
                               mu1 = mixmdl$mu[idx1], mu2 = mixmdl$mu[idx2],
                               sd1 = mixmdl$sigma[idx1], sd2 = mixmdl$sigma[idx2],
                               rho1 = mixmdl$lambda[idx1], rho2 = mixmdl$lambda[idx2]),
                silent = TRUE)

    if (!"try-error" %in% methods::is(root)) {
      # if (verbose) {
      #   abline(v = root$root, col = "red")
      #   abline(v = mixmdl$mu[idx1] + qnorm(0.99) * mixmdl$sigma[idx1], col = "blue")
      # }
      t <- root$root
    }else{
      t <- 0
    }

  }else{
    t <- 0
  }

  return(t)
}


#' @importFrom stats dnorm

funMixModel <- function(x, mu1, mu2, sd1, sd2, rho1, rho2) {

  stats::dnorm(x, mean = mu1, sd = sd1) * rho1 - stats::dnorm(x, mean = mu2, sd = sd2) * rho2


}
