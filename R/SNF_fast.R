#' Fast similarity network fusion method
#' @param Wall List of matrices. Each element of the list is a square, symmetric matrix that shows affinities of the data points from a certain view.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t Number of iterations for the diffusion process.
#' @return W is the overall status matrix derived
#' @useDynLib CiteFuse
#'
#' @import Rcpp RcppEigen
#'
#' @export
#'

SNF_fast <- function(Wall, K = 20, t = 20)
{

  wall.name.check <- check_wall_names(Wall)
  wall.names <- dimnames(Wall[[1]])
  if (!wall.name.check) {
    warning("Dim names not consistent across all matrices in Wall.\n            Returned matrix will have no dim names.")
  }
  LW <- length(Wall)

  newW <- vector("list", LW)
  nextW <- vector("list", LW)
  # for (i in 1:LW) {
  #   Wall[[i]] <- normalize(Wall[[i]])
  #   Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2
  # }

  Wall <- lapply(1:LW, function(i) {
    w <- normalize(Wall[[i]])
    w <- (w + t(w))/2
  })

  # for (i in 1:LW) {
  #   newW[[i]] <- (.dominateset(Wall[[i]], K))
  # }

  newW <- lapply(1:LW, function(i) {
    (.dominateset(Wall[[i]], K))
  })


  for (i in 1:t) {
    for (j in 1:LW) {
      sumWJ <- matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      for (k in 1:LW) {
        # weights may be added here
        if (k != j) {
          sumWJ <- sumWJ + Wall[[k]]
        }
      }
      nextW[[j]] <- eigenMapMatMult(newW[[j]], (sumWJ/(LW - 1)), t(newW[[j]]))
    }

    ## this for loop can be parallel
    # for (j in 1:LW) {
    #   Wall[[j]] <- normalize(nextW[[j]])
    #   Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2
    # }

    Wall <- lapply(1:LW, function(j) {
      w <- normalize(nextW[[j]])
      w <- (w + t(w))/2
    })
  }
  W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
  for (i in 1:LW) {
    W <- W + Wall[[i]]
  }
  W <- W/LW
  W <- normalize(W)
  W <- (W + t(W))/2
  if (wall.name.check) {
    dimnames(W) <- wall.names
  }
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


normalize <- function(X) {
  row.sum.mdiag <- rowSums(X) - diag(X)
  row.sum.mdiag[row.sum.mdiag == 0] <- 1
  X <- X/(2 * (row.sum.mdiag))
  diag(X) <- 0.5
  return(X)
}

check_wall_names <- function(Wall) {
  name_match <- function(names_A, names_B) {
    return(identical(dimnames(names_A), dimnames(names_B)))
  }
  return(all(unlist(lapply(Wall, FUN = name_match, Wall[[1]]))))
}
