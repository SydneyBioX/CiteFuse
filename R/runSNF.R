#' CiteFuse
#'
#' A function to runSNF for CITE seq data
#'
#'
#' @param sce a SingleCellExperiment
#' @param altExp_name expression name of ADT matrix
#' @param W_list affinity list, if it is NULL, the function will calculate it.
#' @param gene_select whether highly variable genes will be selected
#' for RNA-seq to calcualte simlarity matrix using `scran` package
#' @param dist_cal_RNA similarity metrics used for RNA matrix
#' @param dist_cal_ADT similarity metrics used for ADT matrix
#' @param ADT_subset A vector  indicates the subset that will be used.
#' @param K_knn Number of nearest neighbours
#' @param K_knn_Aff Number of nearest neighbors for computing affinity matrix
#' @param sigma Variance for local model for computing affinity matrix
#' @param t Number of iterations for the diffusion process.
#' @param metadata_names A vector indicates the names of metadata returned
#' @param verbose whether print out the process
#' @param FDR false discovery rate threshold for highly variable gene selection
#' @param bio biological component of the variance threshold for highly 
#' variable gene selection 
#' (see `modelGeneVar` in `scran` package for more details)
#'
#' @return A SingleCellExperiment object with fused matrix results stored
#'
#' @examples
#' data("sce_ctcl_subset", package = "CiteFuse")
#' sce_ctcl_subset <- CiteFuse(sce_ctcl_subset)
#'
#' @importFrom SingleCellExperiment SingleCellExperiment logcounts
#' @importFrom Matrix rowSums
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SNFtool affinityMatrix SNF
#' @importFrom propr propr
#' @importFrom stats cor
#' 
#' @references 
#' 
#' B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, 
#' A Goldenberg (2014) Similarity Network Fusion: a fast and effective method 
#' to aggregate multiple data types on a genome wide scale. 
#' Nature Methods. Online. Jan 26, 2014
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
                     K_knn_Aff = 30,
                     sigma = 0.45,
                     t = 20,
                     metadata_names = NULL,
                     verbose = TRUE,
                     FDR = 0.05,
                     bio = 0.01) {


    if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
        stop("sce does not contain logcounts...
               we will perform normalize() to get logcounts")
    }

    if (is.null(W_list)) {
        if (gene_select) {
            hvg <- selectHVG(sce,
                             FDR = FDR,
                             bio = bio)
        } else {
            hvg <- seq_len(nrow(sce))
        }



        rna_mat <- as.matrix(SingleCellExperiment::logcounts(sce)[hvg,])

        if (dist_cal_RNA == "correlation") {
            dist_rna <- as.matrix(1 - stats::cor(rna_mat))
        }

        adt_mat <- assay(altExp(sce, altExp_name), "counts")

        if (is.null(ADT_subset)) {
            ADT_subset <- rownames(adt_mat)
        }

        adt_dist <- suppressMessages(propr(as.matrix(adt_mat[ADT_subset, ])))



        if (dist_cal_ADT == "propr") {
            dist_adt <- 1 - as.matrix(adt_dist@matrix)
        }

        if (verbose) {
            cat("Calculating affinity matrix \n")
        }

        W1 <- SNFtool::affinityMatrix(dist_adt, K = K_knn_Aff, sigma = sigma)
        W2 <- SNFtool::affinityMatrix(dist_rna, K = K_knn_Aff, sigma = sigma)
        W_list <- list(W1, W2)
    }

    if (verbose) {
        cat("Performing SNF  \n")
    }

    W <- SNFtool::SNF(W_list, K = K_knn, t = t)

    if (is.null(metadata_names)) {
        metadata_names <- c("SNF_W", "ADT_W", "RNA_W")
    }

    S4Vectors::metadata(sce)[[metadata_names[1]]] <- W
    S4Vectors::metadata(sce)[[metadata_names[2]]] <- W1
    S4Vectors::metadata(sce)[[metadata_names[3]]] <- W2

    return(sce)
}



#' @importFrom scran modelGeneVar

selectHVG <- function(sce,
                      FDR = 0.05,
                      bio = 0.01) {

    decomp <- scran::modelGeneVar(sce)
    decomp <- decomp[order(decomp$bio, decreasing = TRUE), ]
    hvg <- rownames(decomp)[decomp$FDR < FDR & decomp$bio > bio]

    return(hvg)
}
