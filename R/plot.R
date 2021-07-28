#I hate this but we're gonna make a call to globalVariables here to satisfy R cmd check.
utils::globalVariables(c("from", "."))

#' Plot subnetwork identified using the compute_crosstalk function
#'
#' Convenience function for plotting crosstalkers - if you want to make more
#'     customized/dynamic figures, there are lots of packages that can facilitate that,
#'     including: \code{visnetwork}, \code{ggraph}, and even the base R plotting library
#'
#' @param crosstalk_df a dataframe containing the results of \code{compute_crosstalk}
#' @param label_prop Proportion of nodes to label - based on degree
#' @param prop_keep How many proteins do we want to keep in the visualization (as a proportion of total) -
#'      subsets on top x proteins ranked by affinity score
#' @inheritParams compute_crosstalk
#'
#' @importFrom rlang .data
#'
#' @export

plot_ct <- function(crosstalk_df, g, label_prop = 0.1,
                    prop_keep = 0.4) {

  #make sure g is an igraph object
  if(!igraph::is.igraph(g)) {
    stop("g must be an igraph object.")
  }

  #make sure input is valid compute_crosstalk output
  check_crosstalk(crosstalk_df = crosstalk_df)

  crosstalk_df <- dplyr::slice_max(crosstalk_df,
                                   order_by = .data$affinity_score,
                                   prop = prop_keep)

  #generate vector of seeds
  seeds_df <- dplyr::filter(crosstalk_df, .data$seed == "yes")
  seed_proteins <- seeds_df$gene_id

  #make subgraph of g for a given crosstalk_df
  g_ct <- crosstalk_subgraph(crosstalk_df = crosstalk_df, g = g,
                             seed_proteins = seed_proteins)

  #igraph::tkplot(g_ct)
  ggraph::ggraph(g_ct) +
    ggraph::geom_node_point(ggplot2::aes(size = .data$degree,
                                         color = .data$seed_label)) +
    ggraph::geom_edge_fan(alpha = 0.4, color = "blue") +
    ggraph::geom_node_label(ggplot2::aes(label =
                                           ifelse(.data$degree_rank > 1-label_prop, .data$name, "")),
                            repel = TRUE, hjust = 2, size = 4)

}


#' Check to make sure incoming object is a valid crosstalk df.
#'
#' This function is a helper function for \code{plot_ct} that verifies the input is a valid output of compute_crosstalk
#' @inheritParams plot_ct
#' @inheritParams compute_crosstalk
#'

check_crosstalk <- function(crosstalk_df) {
  #make sure it is a dataframe
  if(!is.data.frame(crosstalk_df)) {
    stop("crosstalk_df must be a valid output of compute_crosstalk")
  }

  #make sure columns match up
  target_cols = c("gene_id", "mean_p", "var_p", "nobs", "seed",
                  "affinity_score", "Z", "p_value", "adj_p_value")

  if(!all(target_cols %in% colnames(crosstalk_df))) {
    stop("column names do not match what is expected")
  }

}

#' Helper function to generate subgraph from crosstalk_df output of \code{compute_crosstalk}
#'
#' Useful if the user wants to carry out further analysis or design custom visualizations.
#'
#' @importFrom magrittr %>%
#'
#' @inheritParams plot_ct
#' @inheritParams compute_crosstalk
#'
#' @return a tidygraph structure containing information about the crosstalkr subgraph
#'
#' @export

crosstalk_subgraph <- function(crosstalk_df, g, seed_proteins) {
  #verify that crosstalk_df is a valid compute_crosstalk output
  check_crosstalk(crosstalk_df = crosstalk_df)

  #if g isn't an igraph object this will fly an error.
  g <- igraph::induced_subgraph(g, v = crosstalk_df$gene_id)

  #we only want to keep edges that attach to a seed protein - is this really true?
  #seed_edges <- igraph::E(g)[ from(seed_proteins)]

  #igraph::subgraph.edges(g, eids = seed_edges) %>%
  g <-  g %>% tidygraph::as_tbl_graph() %>%
    tidygraph::mutate(
      degree = igraph::degree(.),
      degree_rank = dplyr::percent_rank(.data$degree),
      seed_label = ifelse(.data$name %in% seed_proteins, "seed", "crosstalker"))

}


