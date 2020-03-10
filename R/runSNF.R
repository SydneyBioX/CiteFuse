#' CiteFuse
#'
#' A function to runSNF for CITE seq data
#'
#'
#' @param sce a SingleCellExperiment#' @param return_sce A logical input indicates whether a \code{SingleCellExperiment}
#' object will be return
#' @param altExp_name expression name of ADT matrix
#' @param W_list affinity list, if it is NULL, the function will calculate it.
#' @param gene_select whether highly variable genes will be selected for RNA-seq to calcualte simlarity matrix
#' @param dist_cal_RNA similarity metrics used for RNA matrix
#' @param dist_cal_ADT similarity metrics used for ADT matrix
#' @param ADT_subset A vector  indicates the subset that will be used.
#' @param K_knn Number of nearest neighbours
#' @param t Number of iterations for the diffusion process.
#' @param metadata_names A vector indicates the names of metadata returned
#' @param verbose whether print out the process
#'
#' @return A SingleCellExperiment object with fused matrix results stored
#'
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom Matrix rowSums
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SNFtool affinityMatrix
#' @importFrom propr propr
#' @importFrom stats cor
#'
#' @export



CiteFuse <- function(sce,
                   altExp_name = "ADT",
                   W_list = NULL,
                   gene_select = TRUE,
                   dist_cal_RNA = "correlation",
                   dist_cal_ADT = "propr",
                   ADT_subset = NULL,
                   K_knn = 20,
                   t = 20,
                   metadata_names = NULL,
                   verbose = TRUE
                   ) {



  if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
    stop("sce does not contain logcounts... we will perform normalize() to get logcounts")
  }




  if (is.null(W_list)) {
    if (gene_select) {
      hvg <- selectHVG(sce)
    } else {
      hvg <- seq_len(nrow(sce))
    }


    # if (!"logcounts" %in% SummarizedExperiment::assayNames(altExp(sce, altExp_name))) {
    #   warning("ADT does not contain logcounts... we will perform normaliseExprs() to get clr counts")
    #   sce <- normaliseExprs(sce, altExp_name, "clr")
    # }


    rna_mat <- as.matrix(SingleCellExperiment::logcounts(sce)[hvg,])

    if (dist_cal_RNA == "correlation") {
      dist_rna <- as.matrix(1 - stats::cor(rna_mat))
    }

    adt_mat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, altExp_name), "counts")

    if (is.null(ADT_subset)) {
      ADT_subset <- rownames(adt_mat)
    }

    adt_dist <- suppressMessages(propr::propr(as.matrix(adt_mat[ADT_subset, ])))



    if (dist_cal_ADT == "propr") {
      dist_adt <- 1 - as.matrix(adt_dist@matrix)
    }

    if (verbose) {
      cat("Calculating affinity matrix \n")
    }

    W1 <- SNFtool::affinityMatrix(dist_adt, K = K_knn, sigma = 0.5)

    W2 <- SNFtool::affinityMatrix(dist_rna, K = K_knn, sigma = 0.5)



    W_list <- list(W1, W2)
  }

  if (verbose) {
    cat("Performing SNF  \n")
  }

  W <- SNF_fast(W_list, K = K_knn, t = 20)

  if (is.null(metadata_names)) {
    metadata_names <- c("SNF_W", "ADT_W", "RNA_W")
  }

  S4Vectors::metadata(sce)[[metadata_names[1]]] <- W
  S4Vectors::metadata(sce)[[metadata_names[2]]] <- W1
  S4Vectors::metadata(sce)[[metadata_names[3]]]<- W2

  return(sce)
}



#' @importFrom scran trendVar decomposeVar

selectHVG <- function(sce,
                      FDR = 0.05,
                      bio = 0.1) {
  fit <- scran::trendVar(sce,
                         use.spikes = FALSE)

  decomp <- scran::decomposeVar(sce, fit)
  decomp <- decomp[order(decomp$bio, decreasing = TRUE), ]
  # top.hvgs <- order(decomp$bio, decreasing=TRUE)
  hvg <- rownames(decomp)[decomp$FDR < FDR & decomp$bio > bio]

  return(hvg)
}
