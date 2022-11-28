# Functions to facilitate permutation testing for all pipelines######

#' Bootstrap null distribution for RWR
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
#' @inheritParams prep_stringdb
#'
#' @importFrom magrittr %>%
#' @return data frame containing mean/ standard deviation for null distribution
#' @examples
#' \donttest{
#' #g <- prep_biogrid()
#' #bootstrap_null(seed_proteins = c("EGFR", "KRAS"), g= g, ncores = 1, n = 10)
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
  #limit to only seed proteins that are actually in g
  if(is.null(rownames(w))) {
    seed_proteins <-seed_proteins[seed_proteins %in% 1:nrow(w)]
  } else {
    seed_proteins <- seed_proteins[seed_proteins %in% rownames(w)]
  }

  seeds <- match_seeds(g = g, seed_proteins = seed_proteins, n = n)


  if(ncores == 1) {
    null_dist <- list()
    for(i in 1:(n/agg_int)) {
        agg_df <- list()
        for(j in 1:agg_int) {
          counter <- (i-1)*agg_int + j #this keeps us in line with the number of entries in the "seeds" list
          seeds_i <- unlist(seeds[[counter]])


          if(all(!grepl("\\D", seeds_i))) { #coerce the seeds to numbers if possible
            seeds_i <- as.numeric(seeds_i)
          }
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
        null_dist[[i]] <- agg_df
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

          if(all(!grepl("\\D", seeds_i))) { #coerce the seeds to numbers if possible
            seeds_i <- as.numeric(seeds_i)
          }
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
#' @importFrom stats approxfun density
#'
#' @return list of character vectors: randomly generated seed proteins with a similar
#'         degree distribution to parent seed proteins

match_seeds <- function(g, seed_proteins, n, set_seed = NULL) {

  #make sure seed_proteins are a character vector of length >=1
  if(!is.character(seed_proteins[1])) {
    seed_proteins <- as.character(seed_proteins)
  }
  if(is.numeric(set_seed)) {
    oldseed <- .Random.seed
    set.seed(set_seed)
    withr::defer(.Random.seed <- oldseed) #clean up so other functions relying on RNG are not affected.
  }
  num_seeds <- length(seed_proteins)

  #grab degree for the seeds and for the rest of the graph
  degree_seeds <- (degree = igraph::degree(g, v = seed_proteins))
  degree_all <- igraph::degree(g)

  #assign sample seeds vector... use names
  if(is.null(names(degree_all))) {
    names(degree_all) <- as.character(1:length(degree_all))
  }

  #make sure no seeds not in g make it in
  seed_proteins <- seed_proteins[seed_proteins %in% names(degree_all)]

  #now make sure no seed_proteins are in degree_all
  degree_all <- degree_all[!(names(degree_all) %in% seed_proteins)]

  #grab density distribution of seeds
  if(num_seeds != 1) {
    seeds_density <- density(degree_seeds, n = length(degree_all))
    predict_density <- approxfun(seeds_density) #function to generate probabilities based on seed distribution
    prob_dist <- predict_density(degree_all)
    prob_dist[is.na(prob_dist)] <- 0
  } else {
    degree_all <- degree_all[degree_seeds -3 < degree_all & degree_all < degree_seeds + 3] # this line of code limits the sample-able space to +- 3 of the input degree
  }

  #create n seed sets with a similar degree distribution
  sample_seeds <- list()
  for(i in 1:n) {
    if(num_seeds != 1) {
      samp <- sample(degree_all, size = num_seeds, prob = prob_dist)
    } else {
      samp <- sample(degree_all, size = num_seeds) #just take a random seed with a degree +- 3 of the seed
    }
    sample_seeds[[i]] <- names(samp)
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
#' @returns data.frame containing summary statistics for the computed null distribution
#' @export

dist_calc <- function(df, seed_proteins) {

  #pivot longer to prep for summarise
  null_dist <- tidyr::pivot_longer(df, cols = tidyr::everything(), names_to = "node", values_to = "p")

  if(!all(stringr::str_detect(null_dist$node, "\\D"))) {
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
#' @returns data.frame

final_dist_calc <- function(df_list) {
  df <- df_list %>% dplyr::bind_rows() %>%
    dplyr::group_by(.data$node) %>%
    dplyr::summarise(mean_p = mean(.data$mean_p),
                     var_p = mean(.data$var_p),
                     nobs = sum(.data$nobs),
                     seed = unique(seed))
  return(df)
}


#' function to compute null distribution of dnp
#'
#' `compute_null_dnp` calculates a null distribution for the change in network potential for
#'  for each node in a cell signaling network.
#'
#'  The input for this function will be the output of [compute_dnp()].
#'   To compute the null distribution, the nodes in the provided cell signaling
#'   network will be randomly permuted `n` times, with dnp computed or each new
#'   cell signaling network. The mean and standard error of dnp for this set of random
#'   networks will constitute the null model that we will use for comparison.
#'   Be warned that this operation is extremely expensive computationally. It is
#'   recommended to either use a high-performance cluster or limit the computation of the
#'   null distribution to a small number of nodes.
#'   To distribute the workload over multiple cores, just specify ncores.
#'
#' @seealso [compute_dnp()] and [compute_np()]
#'
#' @param df output of [compute_dnp()]
#' @param n number of permutations
#' @param n_genes integer describing number of genes per sample that we will compute the null distribution for
#' @inheritParams compute_dnp
#'
#' @importFrom foreach %dopar% %do% %:%
#' @importFrom magrittr %>%
#'
#' @export
#' @returns df, also saves to cache if specified

compute_null_dnp <- function(cache = NULL, df, ppi = "biogrid", n,
                         n_genes = 50, experiment_name, ncores = 4,
                         min_score = NULL) {
  if(is.null(cache)) {
    stop("please provide a cache for saving and loading data/output")
  }
  #just return the file if we've already done this
  if(file.exists(paste0(cache, experiment_name, "dnpNull.Rda"))) {
    load(paste0(cache, experiment_name, "dnpNull.Rda"))
    return(null_df)
  }

  #load ppi
  g <- load_ppi(cache = cache, min_score = min_score, ppi = ppi)

  samples <- unique(df$sample_name)

  #grab top n genes per sample to iterate over
  v_df = get_topn(df =df, n_genes = n_genes)

  #get rid of dnp
  df <- df %>% dplyr::select(-.data$dnp)
  #do we want to do a grouped split on sample and then iterate directly through the list? probably
  df_list <- df %>% dplyr::group_by(.data$sample_name) %>%
    dplyr::group_split()

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  null_df <-
    foreach::foreach(j = iterators::iter(df_list),
                     .packages = "crosstalkr",
                     .combine = "rbind"
    ) %dopar% {


      mem_df <- list()
      agg_df <- list()
      z <- 1
      for(i in 1:n) {

        #get the info on the top genes for a given sample
        j_sample <- j$sample_name[1]
        v_j <- v_df[,j_sample][[1]][[1]] #Indexing is gross but this just grabs the top n genes for a given sample
        #scramble the graph
        g_i <- get_random_graph(g)

        #calc change in network potential
        dnp_ij <- calc_dnp_i(j, g_i, v_rm = v_j, keep_all = FALSE)

        #for testing only
        # dnp_ij <- dplyr::filter(j, .data$gene_name %in% v_j)
        # dnp_ij$dnp <- rnorm(nrow(dnp_ij), mean = 50, sd = 20)

        mem_df[[i]] <- dnp_ij

        if(i %% 10 == 0) {
          agg_df[[z]] <- combine_null(mem_df)
          agg_df[[z]]$z <- z
          mem_df <- list() #clear "mememory
          z <- z + 1#progress ticker
        }
      }
      return(final_combine(agg_df)) #pass back to the main foreach loop

    }
  # null_j_df <- cell_line_df %>%
  #   dplyr::group_by(.data$gene_name, .data$sample_name) %>%
  #   dplyr::summarise(mean_dnp = mean(.data$np),
  #                    sd_dnp = sd(.data$dnp),
  #                    n = dplyr::n())
  #return(null_j_df)



  #close out parallel execution
  parallel::stopCluster(cl)

  #save
  save(null_df, file = paste0(cache, experiment_name, "dnpNull.Rda"))

  #return
  return(null_df)
}

#' Helper function for compute_null_dnp - returns a graph with randomly permuted edges.
#'
#' currently just a wrapper for igraph::rewire but may add more functionality in the future
#'
#' @param g graph to be permuted
#'
#' @seealso [igraph::rewire()]
#' @export
#' @returns igraph

get_random_graph <- function(g) {
  g2 <- igraph::rewire(g, igraph::keeping_degseq(niter = igraph::vcount(g)*100))
  return(g2)
}

#' Helper function for compute_null_dnp - returns the top n genes by dnp for each sample
#' @inheritParams compute_null_dnp
#' @importFrom magrittr %>%

get_topn <- function(df, n_genes) {
  df %>% dplyr::group_by(.data$sample_name) %>%
    dplyr::slice_max(order_by = .data$dnp, n = n_genes) %>%
    dplyr::select(.data$gene_name, .data$sample_name) %>%
    tidyr::pivot_wider(names_from = .data$sample_name,
                       values_from = .data$gene_name,
                       values_fn = list)
}
#' .combine function for compute_null foreach looping structure
#'
#' @param x aggregated data structure
#' @importFrom magrittr %>%
#' @export
#' @returns data.frame

combine_null <- function(x) {

  #don't want to try and compute sd with one sample... thats where the trouble is ocming I t
  x <- dplyr::bind_rows(x) %>%
    dplyr::group_by(.data$gene_name, .data$sample_name) %>%
    dplyr::summarise(mean_dnp = mean(.data$dnp),
                     sd_dnp = sd(.data$dnp, na.rm = TRUE),
                     n = dplyr::n())

  return(x)
}

#' final .combine function to run in compute_null_dnp foreach looping structure
#'
#' @param x aggregated info
#' @importFrom magrittr %>%
#' @export
#'
#' @returns data.frame
final_combine <- function(x) {

  x <- dplyr::bind_rows(x)
  x <- x %>%
    dplyr::group_by(.data$gene_name, .data$sample_name) %>%
    dplyr::summarise(mean_dnp = mean(.data$mean_dnp),
                     sd_dnp = mean(.data$sd_dnp, na.rm = TRUE),
                     n = sum(.data$n))
  return(x)
}



