#' Perform random walk with repeats on a sparse matrix
#'
#' This function borrows heavily from the RWR function in the RANKS package (cite here)
#'
#' @export
#'
#' @param seed_proteins user defined seed proteins
#' @param w  The adjacency matrix of a given graph in sparse format - dgCMatrix
#' @param gamma restart probability
#' @param eps maximum allowed difference between the computed probabilities at the steady state
#' @param norm if True, w is normalized by dividing each value by the column sum.
#' @param tmax the maximum number of iterations for the RWR
#'
sparseRWR <- function(w, seed_proteins, gamma = 0.6, eps = 1e-10, tmax = 1000, norm = TRUE) {

  # divide the values in each column by the within-column sum
  if(norm == TRUE) {
    w <- norm_colsum(w)
  }

  ##The code below is reproduced from the RWR function of the RANKS package.
  #setup intial p vector
  n <- nrow(w)
  p0 <- p <- numeric(n)
  names(p) <- names(p0) <- rownames(w)
  num_seeds <- length(seed_proteins)
  p0[seed_proteins] <- 1/num_seeds
  p <- p0

  for (t in 1:tmax) {
    # pi+1 = (1 − d)Wpi + dr
    pold <- p
    p <- ((1-gamma) * as.vector(pold %*% w)) + gamma * p0

    #check if the exit condition has been met.
    if (norm1(p-pold) < eps) {
      break
    }
  }
  return(list(p=p, seed_proteins=seed_proteins, n.iter=t));

}

#' Function to normalize adjacency matrix by dividing each value by the colsum.
#'
#' @inheritParams sparseRWR
#'
norm_colsum <- function(w) {
  colsums <- Matrix::colSums(w)
  w <- w/colsums
  return(w)
}


#' Function that computes the norm 1 of a numeric vector
#'
#' This function is reproduced from the source code of the RANKS package (not exported).
#'
#' @param x : numeric vector
#' @return a single real value (the norm1 of the input vector)
norm1 <- function(x) {
  return(sum(abs(x)));
}


