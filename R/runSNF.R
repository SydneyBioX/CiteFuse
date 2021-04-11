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
#' @param topN top highly variable genes are used 
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
                     t = 10,
                     metadata_names = NULL,
                     verbose = TRUE,
                     topN = 2000) {


    if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
        stop("sce does not contain logcounts...
               we will perform normalize() to get logcounts")
    }

    if (is.null(W_list)) {
        if (gene_select) {
            hvg <- selectHVG(sce,
                             topN = topN)
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

        W1 <- affinityMatrix(dist_adt, K = K_knn_Aff, sigma = sigma)
        W2 <- affinityMatrix(dist_rna, K = K_knn_Aff, sigma = sigma)
        W_list <- list(W1, W2)
    }

    if (verbose) {
        cat("Performing SNF  \n")
    }

    W <- SNF(W_list, K = K_knn, t = t)

    if (is.null(metadata_names)) {
        metadata_names <- c("SNF_W", "ADT_W", "RNA_W")
    }

    S4Vectors::metadata(sce)[[metadata_names[1]]] <- W
    S4Vectors::metadata(sce)[[metadata_names[2]]] <- W1
    S4Vectors::metadata(sce)[[metadata_names[3]]] <- W2

    return(sce)
}



#' @importFrom scran modelGeneVar getTopHVGs

selectHVG <- function(sce,
                      topN = 2000) {

    decomp <- scran::modelGeneVar(sce)
    hvg <- getTopHVGs(decomp, n = topN)
    return(hvg)
}



SNF <- function(Wall, K=20, t=20) {
    # Similarity Network Fusion takes multiple views of a network (Wall) and
    # fuses them together to create a overall affinity matrix.
    #
    # Args:
    #   Wall: List of matrices, each element is a square symmetric affinity 
    #       matrix.
    #   K: Number of neighbors used in the K-nearest neighbours step,??? more details???
    #   t: Number of iterations for the diffusion process
    #
    # Returns:  
    #   W: Unified similarity graph of all data types in Wall. 

    
    #Check if Wall names are consistant across all matrices in Wall
    wall.name.check <- check_wall_names(Wall)
    wall.names <- dimnames(Wall[[1]])
    if(!wall.name.check){
        warning("Dim names not consistent across all matrices in Wall.
            Returned matrix will have no dim names.")
    }
    
    LW <- length(Wall)
    

    
    #Normalize different networks to avoid scale problems.
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    for(i in 1:LW){
        Wall[[i]] <- normalize(Wall[[i]])
        Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2
    }
    
    ### Calculate the local transition matrix. (KNN step?)
    for(i in 1:LW){
        newW[[i]] <- (.dominateset(Wall[[i]], K))
    }
    
    #Perform the diffusion for t iterations
    for (i in 1:t) {
        for(j in 1:LW){
            sumWJ <- matrix(0,dim(Wall[[j]])[1], dim(Wall[[j]])[2])
            for(k in 1:LW){
                if(k != j) {
                    sumWJ <- sumWJ + Wall[[k]]
                }
            }
            nextW[[j]] <- newW[[j]] %*% (sumWJ/(LW-1)) %*% t(newW[[j]])
        }
        
        #Normalize each new obtained networks.
        for(j in 1 : LW){
            Wall[[j]] <- normalize(nextW[[j]])
            Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2;
        }
    }
    
    # Construct the combined affinity matrix by summing diffused matrices
    W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
    for(i in 1:LW){
        W <- W + Wall[[i]]
    }
    
    W <- W/LW
    W <- normalize(W)
    W <- (W + t(W)) / 2
    
    if(wall.name.check){
        dimnames(W) <- wall.names
    } 
    
    return(W)  
}


affinityMatrix <- function(diff,K=20,sigma=0.5) {
    # Computes the affinity matrix for a given distance matrix
    # 
    # Args:
    #   diff: Distance matrix 
    #   K: Number of nearest neighbours to sparsify similarity
    #   sigma: Variance for local model
    #
    # Returns:
    #   Affinity matrix using exponential similarity kernel scaled by k nearest
    #       neighbour average similarity
    #
    
    N <- nrow(diff)
    
    diff <- (diff + t(diff)) / 2
    diag(diff) <- 0
    sortedColumns <- as.matrix(t(apply(diff,2,sort)))
    
    finiteMean <- function(x){
        return(mean(x[is.finite(x)]))
    }
    means <- apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
    
    avg <- function(x,y){
        return((x+y)/2)
    }
    Sig <- outer(means,means,avg)/3*2 + diff/3 + .Machine$double.eps;
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(diff, 0, sigma*Sig, log = FALSE)
    
    W <- (densities + t(densities)) / 2
    return(W)
}



.dominateset <- function(xx,KK=20) {
    ###This function outputs the top KK neighbors.	
    
    zero <- function(x) {
        s = sort(x, index.return=TRUE)
        x[s$ix[1:(length(x)-KK)]] = 0
        return(x)
    }
    normalize <- function(X) X / rowSums(X)
    A = matrix(0,nrow(xx),ncol(xx));
    for(i in 1:nrow(xx)){
        A[i,] = zero(xx[i,]);
        
    }
    
    
    return(normalize(A))
}



check_wall_names <- function(Wall){
    # Checks if dimnames are consistant across all matrices in Wall
    #   #Move to internal functions?
    # Args:
    #   Wall: List of matrices
    # Returns:
    #   logical: True/False indicator of dimnames equivalence
    name_match <- function(names_A, names_B){
        return(identical(dimnames(names_A), dimnames(names_B)))
    }
    
    return(all(unlist(lapply(Wall, FUN=name_match, Wall[[1]]))))
}

#Normalization method for affinity matrices
normalize <- function(X){
    row.sum.mdiag <- rowSums(X) - diag(X) 
    #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
    row.sum.mdiag[row.sum.mdiag == 0] <- 1   
    X <- X/(2*(row.sum.mdiag))
    diag(X) <- 0.5
    return(X)
}