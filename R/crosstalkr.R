#' crosstalkr: A package for the identification of functionally relevant subnetworks from high-dimensional omics data.
#'
#' crosstalkr provides a key user function, \code{compute_crosstalk} as well
#'  as several additional functions that assist in setup and visualization (under development).
#'
#' @section crosstalkr functions:
#'
#' \code{compute_crosstalk} calculates affinity scores of all proteins in a protein-protein
#'  interaction network relative to user-provided seed proteins.
#'
#' \code{sparseRWR} performs random walk with restarts on a sparse matrix.
#'  Compared to dense matrix implementations, this should be extremely fast.
#'
#' \code{bootstrap_null} Generates a null distribution based on n calls to
#'  \code{sparseRWR}
#'
#' \code{setup_init} manages download and storage of interactome data to
#'  speed up future analysis
#'
#'
#' @docType package
#'
#' @name crosstalkr
#'
NULL
#>NULL
