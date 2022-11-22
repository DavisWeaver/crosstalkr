### Functions to help with users trying to include these network analyses in a pipeline ######
# These functions will allow users to provide an expression matrix / tidy data frame and then return a tidy data frame with network potential/ change in np

#' main function to compute np from a user-provided expression matrix.
#'
#' @param cache user-provided filepath for where to store data etc
#' @param ppi should we use biogrid or stringdb for the PPI
#' @param mir_paper are we running this in the context of the mir paper? a few quirks of that data
#' @param exp_mat expression matrix where columns are samples and rows are features
#' @param ncores number of cores to use for calculations
#' @param min_score if ppi is stringdb, which mininum score should we use to filter edges?
#' @param experiment_name name of the experiment for saving output.
#' @return tidy data frame with one column for expression and another for np
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom stats sd
#' @importFrom utils download.file unzip
#' @export
#'
#'

compute_np <- function(cache = NULL, experiment_name, ppi = "biogrid", min_score = NULL,
                       exp_mat, mir_paper = TRUE, ncores = 1) {

  if(is.null(cache)) {
    stop("please provide a cache for saving and loading data/output")
  }
  #just return the file if we've already done this
  if(file.exists(paste0(cache, experiment_name, "np.Rda"))) {
    load(paste0(cache, experiment_name, "np.Rda"))
    return(df_np)
  }
  #load ppi
  g <- load_ppi(cache = cache, min_score = min_score, ppi = ppi)

  #convert expression matrix to tidy data frame + do some cleaning
  df <- tidy_expression(exp_mat)

  if(mir_paper == TRUE) { #just split up the information in the sample_id
    df <- experiment_breakout(df)
  }


  samples <- unique(df$sample_name)

  #set up parallel execution
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  out_list <-
    foreach::foreach(i = 1:length(samples)) %dopar% {
      #isolate one cell line
      df_i <- df %>%
        dplyr::filter(.data$sample_name == samples[i])



      #add np to df_i for that cell line
      disruptr::calc_np_i(df_i, g = g)

    }
  #close out parallel execution
  parallel::stopCluster(cl)
  df_np = dplyr::bind_rows(out_list)

  save(df_np, file = paste0(cache, experiment_name, "np.Rda"))

  return(df_np)


}

#' main function to compute delta np for every gene in a given dataframe - assumes compute_np has already been run for a given dataset
#'
#' This function takes a tidy dataframe as input containing RNA sequencing data for one or more samples and conducts in-silico repression.
#' Make sure to run with the same arguments for ppi and cache to maintain consistency for a given pipeline.
#'
#' @inheritParams compute_np
#' @param df dataframe output of compute_np
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @export
#' @returns data.frame

compute_dnp <- function(cache = NULL, df, experiment_name, ppi, ncores = 1,
                        min_score= NULL) {

  if(is.null(cache)) {
    stop("please provide a cache for saving and loading data/output")
  }

  #just return if we have done this before
  if(file.exists(paste0(cache, experiment_name, "dnp.Rda"))) {
    load(paste0(cache, experiment_name, "dnp.Rda"))
    return(df_dnp)
  }
  #load ppi
  g <- load_ppi(cache = cache, min_score = min_score, ppi = ppi)

  #lose any NAs for network potential
  df <- dplyr::filter(df, !is.na(.data$np))

  samples <- unique(df$sample_name)

  #set up parallel execution
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  out_list <-
    foreach::foreach(i = 1:length(samples)) %dopar% {
      #isolate one cell line
      df_i <- df %>%
        dplyr::filter(.data$sample_name == samples[i])

      #add dnp to df_i for that cell line
      disruptr::calc_dnp_i(df_i, g = g)

    }
  #close out parallel execution
  parallel::stopCluster(cl)

  df_dnp = dplyr::bind_rows(out_list)

  save(df_dnp, file = paste0(cache, experiment_name, "dnp.Rda"))

  return(df_dnp)
}

#' helper function to convert expression matrix to tidy dataframe (if not already)
#'
#' @param df dataframe
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @returns data.frame

tidy_expression <- function(df) {

  ##label the first column as gene name
  if(is.factor(unlist(df[,1]))) {
    df[,1] <- as.character(df[,1])
  }
  if(is.character(unlist(df[,1]))) {
    colnames(df)[1] <- "gene_name"
  } else if(is.character(rownames(df))) {
    df$gene_name <- rownames(df)
  }

  #pivot longer
  df <- df %>% tidyr::pivot_longer(-.data$gene_name,
                                   names_to = "sample_name",
                                   values_to = "expression")

  #remove rows with zero expression.
  #don't want log expression here - negative numbers don't play nice
  df <- df %>%
    dplyr::mutate(expression = ifelse(.data$expression < 0, 0, .data$expression))
}

#' helper function to split experiment names into constituent parts
#'
#' this is highly specific  to the miR paper
#'
#' @inheritParams tidy_expression
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @returns data.frame

experiment_breakout <- function(df) {
  experiment_breakout <- stringr::str_split(df$sample_name, pattern = "_",
                                            simplify = TRUE) %>% tibble::as_tibble() %>%
    dplyr::select(1:3)
  colnames(experiment_breakout) <- c("experiment_num", "cell_line", "condition") #breakout the experiment barcode to constituent pieces
  df_ready <- cbind(df, experiment_breakout) %>%
    dplyr::filter(.data$condition == "DMSO") %>%
    dplyr::mutate(gene_name = as.character(.data$gene_name)) #Only working on control cell lines
  return(df_ready)
}

#' helper function to calculate np for one sample
#'
#' @param df dataframe with one cell line + log expression
#' @param g igraph object containing ppi info
#'
#' @return same dataframe with np calculated for each gene.
#'
#' @export

calc_np_i <- function(df, g) {
  #grab expression vector
  exp <- df$expression
  names(exp) <- unlist(df[,1]) #this will break if gene name is not the first column

  #calc np for that expression vector- use cpp internal version
  np <- calc_np_all(exp, g)

  #gene_name must be how we label this columns
  np_df <- data.frame(gene_name = names(np),
                      np = np)

  df <- dplyr::left_join(df, np_df)
  return(df)
}

#' helper function to calculate dnp for one sample
#'
#' @param df dataframe with one cell line + log expression
#' @param g igraph object containing ppi info
#' @param v_rm passed to [node_repression()]
#' @param keep_all logical flag denoting if we should keep genes that we didn't calculate dnp for
#'
#' @return same dataframe with dnp calculated for each gene.
#'
#' @export

calc_dnp_i <- function(df, g, v_rm = NULL, keep_all = TRUE) {
  exp <- df$expression
  names(exp) <- unlist(df[,1])

  if(is.null(v_rm)) {
    v_rm = names(exp)
  }
  dnp_mat <- node_repression(g = g, v_rm = v_rm, exp= exp)

  #sum everything
  dnp <- Matrix::colSums(dnp_mat)

  #add to dataframe and then join with df
  dnp_df <- data.frame(gene_name = names(dnp), dnp = dnp)

  if(keep_all == FALSE) {
    df <- dplyr::filter(df, .data$gene_name %in% v_rm)
  }

  df <- dplyr::left_join(df, dnp_df, by = "gene_name")
  df <- dplyr::mutate(df, dnp = ifelse(is.na(.data$dnp), 0, .data$dnp))

  return(df)
}

#bad: gonna define a bunch of globals to avoid the dreaded R cmd check notes

globalVariables(c("i", "j"))
