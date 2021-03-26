#' Bootstrap null distribution for significance testing
#'
#' This function will generate a bootstrapped null distribution to identify
#' signficant vertices in a PPI given a set of user-defined seed proteins.
#' Bootstrapping is done by performing random walk with repeats repeatedly over "random"
#' sets of seed proteins. Degree distribution of user-provided seeds is used to inform sampling.
#'
#' @importFrom foreach %dopar%
#'
#'
#' @param g igraph object
#' @param n number of random walks with repeats to create null distribution
#' @param set_seed integer to set random number seed - for reproducibility
#' @param ncores Number of cores to use - defaults to 1. Significant speedup can be achieved by using multiple cores for computation.
#' @param seed_name Name to give the cached null distribution - must be a character string
#'
#' @inheritParams sparseRWR
#' @inheritParams setup_init
#'
#' @export

bootstrap_null <- function(seed_proteins, g, n = 1000,
                           gamma=0.6, eps = 1e-10, tmax = 1000,
                           norm = TRUE, set_seed,
                           cache = NULL, seed_name = NULL,
                           ncores = 1) {

  #If file was cached from a previous run - use that- else go through the whole calculation
  if(!is.null(cache)) {
    if(file.exists(paste0(cache, "/", seed_name, "null_dist.Rda"))) {
      load(file = paste0(cache, "/", seed_name, "null_dist.Rda"))
      return(out)
    }
  }
  w <- igraph::as_adjacency_matrix(g) #sparse adjacency matrix.
  #normalize w once now so we don't have to do it repeatedly.
  w <- norm_colsum(w)

  #generate list of degree-similar seed protein vectors.
  seeds <- match_seeds(g = g, seed_proteins = seed_proteins, n = n)

  if(ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    null_dist <-
      foreach::foreach(i = 1:n, .errorhandling = 'pass', .packages = "Matrix") %dopar% {
        seeds_i <- unlist(seeds[[i]])
        crosstalkr::sparseRWR(seed_proteins = seeds_i, w = w, norm = FALSE)[[1]] #norm=FALSE because we already did it.
      }
    parallel::stopCluster(cl)
  } else {
    null_dist <- list()
    for(i in 1:n) {
      seeds_i <- unlist(seeds[[i]])
      null_dist[[i]] <- sparseRWR(seed_proteins = seeds_i, w = w, norm = FALSE)[[1]]
    }
  }

  null_dist <- dplyr::bind_rows(null_dist)
  null_dist <- dist_calc(null_dist,
                         seed_proteins = seed_proteins)

  out <- list(null_dist, seed_proteins)

  if(is.character(cache) & is.character(seed_name)) {
    save(out, file = paste0(cache, "/", seed_name, "null_dist.Rda"))
  } else {
    message("Please provide character string designating a filepath and seed name if you would like to save these results to speed up future analysis.")
  }

  return(out)
}

#' Identify random sets of seeds with similar degree distribution to parent seed proteins
#'
#' This function will generate n character vectors of seeds to be passed to sparseRWR
#' as part of the construction of a boostrapped null distribution for significance testing.
#' @param g igraph object representing the network under study. specified by "ppi" in bootstrap_null
#'
#' @inheritParams bootstrap_null
#'

match_seeds <- function(g, seed_proteins, n, set_seed = NULL) {

  if(is.numeric(set_seed)) {
    oldseed <- .Random.seed
    set.seed(set_seed)
    withr::defer(.Random.seed <- oldseed) #clean up so other functions relying on RNG are not affected.
  }
  num_seeds <- length(seed_proteins)

  #Add degree as a character for each seed_protein.
  degree_seeds <- (degree = igraph::degree(g, v = seed_proteins))


  #make largest possible number of bins relating to the degree_seeds vector
  for(i in seq_along(degree_seeds) - 1) {
    degree_bins <- try(ggplot2::cut_number(degree_seeds, num_seeds - i),
                       silent = TRUE)
    if(!inherits(degree_bins, 'try-error')) {
      break
    }
  }


  #Identify the single number breaks from degree_bins
  bins <- levels(degree_bins)

  #extract the first side of each bin to pass to cut
  breaks <- base::as.numeric(stringr::str_extract(bins, "(?<=\\().*(?=\\,)|(?<=\\[).*(?=\\,)"))
  degree_all <- igraph::degree(g)

  degree_all <- tibble::tibble(gene = names(degree_all),
                               degree = degree_all,
                               degree_bins = cut(degree_all, breaks = breaks))


  #group by breaks
  degree_grouped <- dplyr::group_by(degree_all, degree_bins)

  sample_seeds <- list()

  for(i in 1:n) {
    samp <- dplyr::slice_sample(degree_grouped, n = 1) #n = 1 because the number of groups will always be equal to the number of seeds
    sample_seeds[[i]] <- samp$gene
  }

  return(sample_seeds)

}

#' Internal function that computes the mean/stdev for each gene from a wide-format data frame.
#'
#' This function is called by the high-level function "bootstrap_null".
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats var
#'
#' @param df : numeric vector
#' @inheritParams bootstrap_null
#' @return a 3-column dataframe (gene, )

dist_calc <- function(df, seed_proteins) {

  #pivot longer to prep for summarise
  null_dist <- tidyr::pivot_longer(df, cols = tidyr::everything(), names_to = "gene_id", values_to = "p")

  null_dist <- null_dist %>%
    dplyr::group_by(.data$gene_id) %>%
    dplyr::summarise(mean_p = mean(.data$p),
                     var_p = var(.data$p),
                     nobs = dplyr::n()) %>%
    dplyr::mutate(seed = ifelse(.data$gene_id %in% seed_proteins, "yes", "no"))
  return(null_dist)

}
