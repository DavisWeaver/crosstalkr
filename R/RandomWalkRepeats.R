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
#' @examples
#'# 1) Run Random walk with restarts on a simple matrix
#'v1 = (c(1,1,1,0))
#'v2 = c(0,0,0,1)
#'v3 = c(1,1,1,0)
#'v4 = c(0,0,0,1)
#'w = matrix(data = c(v1,v2,v3,v4), ncol = 4, nrow = 4)
#'sparseRWR(seed_proteins = c(1,3), w = w, norm = TRUE)
#'
#'# 2) Works just as well on a sparse matrix
#'v1 = (c(1,1,1,0))
#'v2 = c(0,0,0,1)
#'v3 = c(1,1,1,0)
#'v4 = c(0,0,0,1)
#'w = matrix(data = c(v1,v2,v3,v4), ncol = 4, nrow = 4)
#'w = Matrix::Matrix(w, sparse = TRUE)
#'sparseRWR(seed_proteins = c(1,4), w = w, norm = TRUE)
#'\donttest{
#'#3) Sample workflow for use with human protein-protein interaction network
#'#g <- prep_biogrid()
#'#w <- igraph::as_adjacency_matrix(g)
#'#sparseRWR(seed_proteins = c("EGFR", "KRAS"), w = w, norm = TRUE)
#'}
#' @return numeric vector, affinity scores for all nodes in graph relative to provided seeds
#'

sparseRWR <- function(seed_proteins, w, gamma = 0.6, eps = 1e-10, tmax = 1000,
                      norm = TRUE) {

  #limit to only seed proteins that are actually in g
  if(is.null(rownames(w))) {
    seed_proteins <-seed_proteins[seed_proteins %in% 1:nrow(w)]
  } else {
    seed_proteins <- seed_proteins[seed_proteins %in% rownames(w)]
  }

  #coerce to a sparse matrix
  w <- Matrix::Matrix(w, sparse = TRUE)

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
  w = Matrix::t(w) #moving this out of the loop so we don't have to do it multiple times

  #Automated check to see if n is the same length as p


  for (t in 1:tmax) {
    # pi+1 = (1 − d)Wpi + dr
    pold <- p

    if(n != length(p)) {
      stop("p not the same length as w, input graph may be too sparse")
    }

    #check if the exit condition has been met
    p <- ((1-gamma) * as.vector(w %*% pold)) + gamma * p0

    if (sum(abs(p-pold)) < eps) {
      break
    }
  }
  return(list(p=p, seed_proteins=seed_proteins, n.iter=t));

}

#' Function to normalize adjacency matrix by dividing each value by the colsum.
#'
#' @inheritParams sparseRWR
#'
#' @return input matrix, normalized by column sums
#' @examples
#'# 1) Normalize by column sum on a simple matrix
#'v1 = (c(1,1,1,0))
#'v2 = c(0,0,0,1)
#'v3 = c(1,1,1,0)
#'v4 = c(0,0,0,1)
#'w = matrix(data = c(v1,v2,v3,v4), ncol = 4, nrow = 4)
#'norm_colsum(w)
#'
#' @export

norm_colsum <- function(w) {
  sums <- Matrix::colSums(w)
  zeros <- sums[sums==0]

  #it was crashing big time when I had zeros on the wrong side of the divider.
  if(length(zeros) > 0) {
    w <- w[!(rownames(w) %in% names(zeros)), !(colnames(w) %in% names(zeros))]
    sums <- Matrix::colSums(w)
  }

  w <- Matrix::t(Matrix::t(w)/ sums) #tried splitting this into 2 lines so it didn't crash my shit.

  return(w)
}


