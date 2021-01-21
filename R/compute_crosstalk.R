#' Identify proteins with a statistically significant relationship to user-provided seeds.
#'
#' \code{compute_crosstalk} returns a dataframe of proteins that are significantly
#' associated with user-defined seed proteins. These identified "crosstalkers"
#' can be combined with the user-defined seed proteins to identify functionally
#' relevant subnetworks. Affinity scores for every protein in the network are
#' calculated using a random-walk with repeats (\code{sparseRWR}). Significance is
#' determined by comparing these affinity scores to a bootstrapped null distribution
#' (see \code{bootstrap_null}).
#'
#' @param significance_level user-defined signficance level for hypothesis testing
#' @param p_adjust adjustment method to correct for multiple hypothesis testing:
#'     defaults to "bonferroni". see \code{\link{stats::p.adjust.methods}} for other potential
#'     adjustment methods.
#'
#' @inheritParams bootstrap_null
#'
#' @inheritParams sparseRWR
#'
#' @inheritParams setup_init
#'
#'
#' @export

compute_crosstalk <- function(seed_proteins, ppi = "stringdb", n = 1000,
                              gamma=0.6, eps = 1e-10, tmax = 1000,
                              norm = TRUE, set_seed,
                              cache, seed_name = NULL,
                              ncores = 1, significance_level = 0.95,
                              p_adjust = "bonferroni")  {
  if(ppi == "biogrid") {
    g <- prep_biogrid(cache = cache)
  } else if (ppi == "stringdb") {
    g <- prep_stringdb(cache = cache)
  } else {
    stop("ppi must be either 'biogrid' or 'stringdb'")
  }

  w <- igraph::as_adjacency_matrix(g) #sparse adjacency matrix.

  #Compute p given seed proteins
  p_seed <- sparseRWR(w = w, seed_proteins = seed_proteins, gamma = gamma,
                      eps = eps, tmax = tmax, norm = norm)
  p_vec <- p_seed[[1]]
  p_df <- tibble::as_tibble(p_vec, rownames = "gene_id")
  colnames(p_df)[colnames(p_df) == "value"] <- "p_test"

  #compute null distribution
  null_dist <- bootstrap_null(seed_proteins = seed_proteins, g = g, n = n,
                              gamma = gamma, eps = eps, tmax = tmax,
                              norm = norm, set_seed = set_seed, cache = cache,
                              seed_name = seed_name, ncores = ncores)
  null_df <- null_dist[[1]]

  df <- dplyr::left_join(null_df, p_df)

  #compute the Z-score and p-value
  df <- dplyr::mutate(df,
                      Z = (p_test - mean_p)/ stdev_p,
                      p_value = 2*pnorm(-abs(Z)),
                      adj_p_value = p.adjust(p_value, method = p_adjust))
  df <- dplyr::filter(df, adj_p_value < 1-significance_level)

  return(df)
}

