#' Bootstrap null distribution for significance testing
#'
#' This function will generate a bootstrapped null distribution to identify
#' signficant vertices in a PPI given a set of user-defined seed proteins.
#' Bootstrapping is done by performing random walk with repeats repeatedly over "random"
#' sets of seed proteins. Degree distribution of user-provided seeds is used to inform sampling.
#'
#' @importFrom foreach %dopar% %do%
#'
#'
#' @param g igraph object
#' @param n number of random walks with repeats to create null distribution
#' @param set_seed integer to set random number seed - for reproducibility
#' @param ncores Number of cores to use - defaults to 1. Significant speedup can be achieved by using multiple cores for computation.
#' @param seed_name Name to give the cached ngull distribution - must be a character string
#' @param agg_int number of runs before we need to aggregate the results - necessary to save memory. set at lower numbers to save even more memory.
#'
#' @inheritParams sparseRWR
#' @inheritParams setup_init
#'
#' @importFrom magrittr %>%
#' @return data frame containing mean/ standard deviation for null distribution
#' @examples
#' \donttest{
#' g <- prep_biogrid()
#' bootstrap_null(seed_proteins = c("EGFR", "KRAS"), g= g, ncores = 1, n = 10)
#' }
#' @export
#'

bootstrap_null <- function(seed_proteins, g, n = 1000, agg_int = 100,
                           gamma=0.6, eps = 1e-10, tmax = 1000,
                           norm = TRUE, set_seed = NULL,
                           cache = NULL, seed_name = NULL,
                           ncores = 1) {

  #agg_int must be <= n
  if(agg_int > n) {
    agg_int <- n
  }
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


  if(ncores == 1) {
    null_dist <-
      foreach::foreach(i = 1:(n/agg_int), .errorhandling = 'pass',
                       .packages = c("Matrix"),
                       .export = c("sparseRWR", "dist_calc",
                                   "norm_colsum")) %do%
      {
        agg_df <- list()
        for(j in 1:agg_int) {
          counter <- (i-1)*agg_int + j #this keeps us in line with the number of entries in the "seeds" list
          seeds_i <- unlist(seeds[[counter]])
          p_i <- sparseRWR(seed_proteins = seeds_i, w = w, norm = FALSE)[[1]]
          if(is.null(names(p_i))) {
            names(p_i) <- as.character(1:length(p_i))
          }
          #norm=FALSE because we already did it.
          agg_df[[j]] <- p_i
          if(j == agg_int) {
            agg_df <- dplyr::bind_rows(agg_df)
            agg_df <- dist_calc(agg_df, seed_proteins = seed_proteins)
            agg_df$run <- i
          }

        }
        return(agg_df)
      }
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    null_dist <-
      foreach::foreach(i = 1:(n/agg_int), .errorhandling = 'pass',
                       .packages = c("Matrix", "magrittr"),
                       .export = c("sparseRWR", "dist_calc",
                                   "norm_colsum")) %dopar%
      {
        agg_df <- list()
        for(j in 1:agg_int) {
          counter <- (i-1)*agg_int + j #this keeps us in line with the number of entries in the "seeds" list
          seeds_i <- unlist(seeds[[counter]])
          p_i <- sparseRWR(seed_proteins = seeds_i, w = w, norm = FALSE)[[1]]
          if(is.null(names(p_i))) {
            names(p_i) <- as.character(1:length(p_i))
          }
          #norm=FALSE because we already did it.
          agg_df[[j]] <- p_i
          if(j == agg_int) {
            agg_df <- dplyr::bind_rows(agg_df)
            agg_df <- dist_calc(agg_df, seed_proteins = seed_proteins)
            agg_df$run <- i
          }

        }
        return(agg_df)
      }
    parallel::stopCluster(cl)
  }


  if(length(null_dist) == 1) {
    null_dist <- null_dist[[1]]
  } else {
    null_dist <- final_dist_calc(null_dist)
  }
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
#' @return list of character vectors: randomly generated seed proteins with a similar
#'         degree distribution to parent seed proteins

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

  #need some kind of dependency where if its not a certain number of nodes then you just randomly select seeds
  if(length(breaks > 1)) {

    degree_df <- tibble::tibble(
      degree = degree_all,
      degree_bins =  cut(degree_all, breaks = breaks))

  } else {
    if(num_seeds == 1) {
      degree_df <- tibble::tibble(
        degree = degree_all,
        degree_bins =  cut(degree_all, breaks = 2))
    } else {
      #if there is only one break - you end up with cut selecting equidistant breaks w/ the value you give for break determing the number.
      degree_df <- tibble::tibble(
        degree = degree_all,
        degree_bins =  cut(degree_all, breaks = num_seeds))
    }

  }

  if(is.null(names(igraph::degree(g)))) {
    degree_df$node <- 1:nrow(degree_df)

  } else {
    degree_df$node <- names(igraph::degree(g))
  }


  #group by breaks
  degree_grouped <- dplyr::group_by(degree_df, degree_bins)

  sample_seeds <- list()

  for(i in 1:n) {
    if(igraph::vcount(g) < 50) {
      sample_seeds[[i]] <- sample(degree_grouped$node, size = num_seeds)
      #so we want to make sure that similar degree distribution is chosen even if the number of seeds is small (lots of errors can occur in that situation)
    } else if (num_seeds == 1 | length(levels(degree_df$degree_bins)) != num_seeds) {
      sample_vec <- degree_df %>%
        dplyr::filter(.data$degree %in% degree_all[seed_proteins]) %>%
        dplyr::select(.data$node) %>% unlist()
      sample_seeds[[i]] <- sample(sample_vec, size = num_seeds)
    } else {
      samp <- dplyr::slice_sample(degree_grouped, n = 1) #n = 1 because the number of groups will always be equal to the number of seeds
      sample_seeds[[i]] <- samp$node
    }
  }

  return(sample_seeds)

}

#' Internal function that computes the mean/stdev for each gene from a wide-format data frame.
#'
#' This function is called by the high-level function "bootstrap_null".
#' Not expected to be used by end-users - we only export it so that environments
#' inside foreach loops can find it.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats var
#'
#' @param df : numeric vector
#' @inheritParams bootstrap_null
#' @return a data frame containing summary statistics for the computed null distribution
#' @export

dist_calc <- function(df, seed_proteins) {

  #pivot longer to prep for summarise
  null_dist <- tidyr::pivot_longer(df, cols = tidyr::everything(), names_to = "node", values_to = "p")

  if(all(stringr::str_detect(null_dist$node, "\\d"))) {
    null_dist$node <- as.numeric(null_dist$node)
  }

  null_dist <- null_dist %>%
    dplyr::group_by(.data$node) %>%
    dplyr::summarise(mean_p = mean(.data$p),
                     var_p = var(.data$p),
                     nobs = dplyr::n()) %>%
    dplyr::mutate(seed = ifelse(.data$node %in% seed_proteins, "yes", "no"))
  return(null_dist)

}

#' Internal function that computes the mean/stdev for each gene from a wide-format data frame.
#'
#' This function is called by the high-level function "bootstrap_null".
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats var
#'
#' @param df_list : list of dataframes from foreach loop in bootstrap_null
#' @return a dataframe

final_dist_calc <- function(df_list) {
  df <- df_list %>% dplyr::bind_rows() %>%
    dplyr::group_by(.data$node) %>%
    dplyr::summarise(mean_p = mean(.data$mean_p),
                     var_p = mean(.data$var_p),
                     nobs = sum(.data$nobs),
                     seed = unique(seed))
  return(df)
}

