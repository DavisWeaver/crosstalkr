#' crosstalkr: A package for the identification of functionally relevant subnetworks from high-dimensional omics data.
#'
#' crosstalkr provides a key user function, \code{compute_crosstalk} as well
#'      as several additional functions that assist in setup and visualization (under development).
#'
#' @section crosstalkr functions:
#'
#' \code{compute_crosstalk} calculates affinity scores of all proteins in a network
#'      relative to user-provided seed proteins. Users can use the human interactome or
#'      provide a network represented as an igraph object.
#'
#' \code{sparseRWR} performs random walk with restarts on a sparse matrix.
#'      Compared to dense matrix implementations, this should be extremely fast.
#'
#' \code{bootstrap_null} Generates a null distribution based on n calls to
#'      \code{sparseRWR}
#'
#' \code{setup_init} manages download and storage of interactome data to
#'      speed up future analysis
#'
#' \code{plot_ct} allows users to visualize the subnetwork identified in
#'      \code{compute_crosstalk}. This function relies on the ggraph framework.
#'      Users are encouraged to use ggraph or other network visualization packages for
#'      more customized figures.
#'
#' \code{crosstalk_subgraph} converts the output of \code{compute_crosstalk} to a
#'      tidygraph object containing only the identified nodes and their connections to the
#'      user-provided seed_proteins. This function also adds degree, degree_rank, and
#'      seed_label as attributes to the identified subgraph to assist in plotting.
#'
#'
#' @docType package
#'
#' @name crosstalkr
#' @
#'
NULL
#>NULL
