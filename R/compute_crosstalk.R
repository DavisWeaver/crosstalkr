#' Identify proteins with a statistically significant relationship to user-provided seeds.
#'
#' \code{compute_crosstalk} returns a dataframe of proteins that are significantly
#'     associated with user-defined seed proteins. These identified "crosstalkers"
#'     can be combined with the user-defined seed proteins to identify functionally
#'     relevant subnetworks. Affinity scores for every protein in the network are
#'     calculated using a random-walk with repeats (\code{sparseRWR}). Significance is
#'     determined by comparing these affinity scores to a bootstrapped null distribution
#'     (see \code{bootstrap_null}). If using non-human PPI from string, refer to the stringdb documentation
#'     for how to specify proteins
#'
#' @param significance_level user-defined signficance level for hypothesis testing
#' @param p_adjust adjustment method to correct for multiple hypothesis testing:
#'     defaults to "holm". see \code{p.adjust.methods} for other potential
#'     adjustment methods.
#' @param use_ppi bool, should g be a protein-protein interaction network? If
#'     false, user must provide an igraph object in \code{g}
#' @param ppi character string describing the ppi to use: currently only "stringdb" and "biogrid" are supported.
#' @param species character string describing the species of interest.
#'     For a list of supported species, see \code{supported_species}.
#'     Non human species are only compatible with "stringdb"
#' @param g igraph network object.
#' @param union bool, should we take the union of string db and biogrid to compute the PPI? Only applicable for the human PPI
#' @param intersection bool, should we take the intersection of string db and biogrid to compute the PPI? Only applicable for the human PPI
#'
#' @inheritParams bootstrap_null
#'
#' @inheritParams sparseRWR
#'
#' @inheritParams prep_stringdb
#'
#' @importFrom stats p.adjust pnorm sd
#' @importFrom utils download.file read.delim unzip
#' @importFrom rlang .data
#'
#' @return data frame containing affinity score, p-value, for all "crosstalkers"
#'         related to a given set of seeds
#'
#' @examples
#' \donttest{
#' #1) easy to use for querying biological networks - n = 10000 is more appropriate for actual analyses
#' #compute_crosstalk(c("EGFR", "KRAS"), n =10)
#' }
#' #2) Also works for any other kind of graph- just specify g (must be igraph formatted as of now)
#' g <- igraph::sample_gnp(n = 1000, p = 10/1000)
#' compute_crosstalk(c(1,3,5,8,10), g = g, use_ppi = FALSE, n = 100)
#'
#'
#' @export

compute_crosstalk <- function(seed_proteins, g = NULL, use_ppi = TRUE,
                              ppi = "stringdb", species = "homo sapiens",
                              n = 1000, union = FALSE, intersection = FALSE,
                              gamma=0.6, eps = 1e-10, tmax = 1000,
                              norm = TRUE, set_seed,
                              cache = NULL, min_score = 700, seed_name = NULL,
                              ncores = 1, significance_level = 0.95,
                              p_adjust = "bonferroni",
                              agg_int = 100)  {

  #check inputs
  if(use_ppi == TRUE){
    if(union & (tolower(species) == "homo sapiens" | as.character(species) == "9606")) {
      g <- ppi_union(cache = cache, min_score = min_score)
    } else if(intersection & (tolower(species) == "homo sapiens" | as.character(species) == "9606")) {
      g <- ppi_intersection(cache = cache, min_score = min_score)
    } else if(ppi == "biogrid" & (tolower(species) == "homo sapiens" | as.character(species) == "9606")) { #first 3 options are only feasible if the species is human
      g <- prep_biogrid(cache = cache)
    } else if (ppi == "stringdb") {
      g <- prep_stringdb(cache = cache, min_score = min_score, species = species)
    } else {
      stop("ppi must be either 'biogrid' or 'stringdb'")
    }
  } else {
    if(!igraph::is.igraph(g)){
      stop("g must be an igraph object")
    }
  }

  w <- igraph::as_adjacency_matrix(g) #sparse adjacency matrix.

  #Compute p given seed proteins
  p_seed <- sparseRWR(seed_proteins = seed_proteins, w = w, gamma = gamma,
                      eps = eps, tmax = tmax, norm = norm)
  p_vec <- p_seed[[1]]
  p_df <- tibble::as_tibble(p_vec, rownames = "node")
  if(!all(stringr::str_detect(p_df$node, "\\D"))) { #if the node names are numeric lets make it so
    p_df$node <- as.numeric(p_df$node)
  }
  colnames(p_df)[colnames(p_df) == "value"] <- "affinity_score"

  #compute null distribution
  null_dist <- bootstrap_null(seed_proteins = seed_proteins, g = g, n = n,
                              gamma = gamma, eps = eps, tmax = tmax,
                              norm = norm, set_seed = set_seed, cache = cache,
                              seed_name = seed_name, ncores = ncores,
                              agg_int = agg_int)
  null_df <- null_dist[[1]]

  df <- dplyr::left_join(null_df, p_df)

  #compute the Z-score and p-value
  df <- dplyr::mutate(df,
                      Z = (.data$affinity_score - .data$mean_p)/ sqrt(.data$var_p),
                      p_value = 2*pnorm(-abs(.data$Z)),
                      adj_p_value = p.adjust(.data$p_value, method = p_adjust))
  df <- dplyr::filter(df, .data$adj_p_value < 1-significance_level)

  return(df)
}

