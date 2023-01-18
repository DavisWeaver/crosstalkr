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
#' @param return_g bool, should we return the graph used? mostly for internal use
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
                              agg_int = 100, return_g = FALSE,
                              network_type = "full")  {

  #check inputs
  if(use_ppi == TRUE){
    g <- load_ppi(ppi = ppi, species = species, min_score = min_score,
                  union = union, intersection = intersection, cache = cache,
                  network_type = network_type)
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

  if(return_g) {
    out <- list(df, g)
  }
  else {
    out <- df
  }
  return(out)
}


#' Generic function to filter either an igraph object or a PPI network
#'
#' @param method str
#' @param g igraph object
#' @param val named numeric vector - some measure of node state (i.e. gene expression in the case of a PPI)
#' @param use_ppi bool - should we use a ppi from online repository?
#' @param igraph_method bool - is the user-provided method an igraph node scoring function?
#' @param n int - number of nodes to include in the returned subgraph
#' @param desc bool - do we want the top or bottom examples of the provided metric
#' @param ... additional params passed to [load_ppi()] or [compute_crosstalk()]
#' @export
#' @returns igraph
#' @seealso [gfilter.ct], [gfilter.np], [gfilter.igraph_method]
#'
gfilter <- function(method=NULL, g = NULL, val = NULL, use_ppi, igraph_method = NULL, n=100, desc=TRUE, ...) {
  select_method <- function(method, igraph_method) {
    if(!is.null(igraph_method)) {
      return(gfilter.igraph_method)
    } else {
      if (method %in% c("ct", "crosstalk", "CT", "compute_crosstalk", "rwr",
                        "RWR")) {
        return(gfilter.ct)
      }
      else if (method %in% c("np", "NP", "network_potential")) {
        return(gfilter.np)
      }
      else if (method %in% c("value", "val", "exp", "expression")) {
        return(gfilter.value)
      }
    }
  }
  func = select_method(method=method, igraph_method=igraph_method)
  if(!is.null(igraph_method)) {
    method = igraph_method
    g <- func(use_ppi = use_ppi, g=g, method=method, n=n, desc=desc, ...)
  } else {
    g <- func(val = val, use_ppi=use_ppi, g=g, n=n, desc=desc, ...)
  }

  return(g)

}
#' Method to filter the graph based on parameters passed to compute_crosstalk
#'
#' @param seeds vector (str or numeric) user provided vertex ids to use as seeds in the random walk with restarts'
#' @param return_df bool should we return a list containing the filtered graph + the RWR output that was used to do the filtering?
#' @param ... additional arguments passed to [compute_crosstalk()]
#'
#' @return igraph object
#' @export
gfilter.ct <- function(seeds, return_df = FALSE, ...) {
  out <- compute_crosstalk(seed_proteins = seeds, return_g = TRUE, ...)
  df <- out[[1]]
  g <- crosstalk_subgraph(df, g=out[[2]], seed_proteins = seeds,tg = FALSE)

  if(return_df) {
    out <- list(g, df)
  } else {
    out <- g
  }
  return(out)
}

#' Method to filter graph based on user provided value
#' @inheritParams gfilter
#' @param val_name str
#' @returns igraph
#'
#' @export

gfilter.value <- function(g, val, use_ppi = TRUE, n = 500, val_name = "value", desc, ...) {
  if(use_ppi) {
    g <- load_ppi(...)
  }

  if(is.null(names(val))) {
    names(val) <- 1:length(val)
  }
  #only include nodes that are in both val and g
  if(is.null(names(igraph::V(g)))) {
    val <- val[names(val) %in% 1:length(igraph::V(g))]
  } else {
    val <- val[names(val) %in% names(igraph::V(g))]
  }

  val = sort(val, decreasing=desc)
  if(n > length(val)) {
    val_keep = val
  } else {
    val_keep = val[1:n]
  }

  #This function actually does the filtering
  g <- add_value(val = val_keep, val_name = val_name, g = g)

  return(g)

}

#' Method to filter graph based on network potential values.
#'
#' convenience function - it just calls gfilter.value after computing np
#'
#' For more information on network potential, see \href{https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1008755}{related paper}
#' @inheritParams gfilter
#' @returns igraph
#' @export
gfilter.np <- function(g, val, use_ppi = TRUE, n = 500, desc, ...) {
  if(use_ppi) {
    g <- load_ppi(...)
  }
  np <- abs(calc_np_all(g=g, exp=val))
  g <- gfilter.value(g=g, val = np, val_name = "np", n=n, desc=desc, use_ppi = FALSE) #we've already loaded the ppi if use_ppi = TRUE
  return(g)
}

#' Method to filter graph based on any igraph method that scores verticies.
#'
#' @inheritParams gfilter
#' @inheritParams gfilter.value
#' @param ... additional parameters passed to [load_ppi]
#' @returns igraph
#' @export


gfilter.igraph_method <- function(g, use_ppi = TRUE, method, n = 500, desc, val_name, ...) {
  if(use_ppi) {
    g <- load_ppi(...)
  }
  if(!is.function(method)) {
    try(attachNamespace("igraph"), silent=TRUE)
    method = try(get(method), silent = TRUE)
    if(inherits(method, 'try-error')) {
      stop("Invalid method specification, please provide a valid igraph method for node scoring")
    }
  }
  val <- method(g)
  g <- gfilter.value(g=g, val=val, val_name = val_name, n=n, use_ppi = FALSE, desc = desc) #we've already loaded the ppi if use_ppi = TRUE
  return(g)
}


