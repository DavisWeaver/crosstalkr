#' function to calculate the network potential for each protein in a user-provided vector
#'
#' Mostly just used to help debug the CPP version - not exported
#' @param exp expression vector - assumed to be a named vector where the values are expression and the names are the gene name
#' @param g igraph object - will be filtered so that only nodes found in both exp and g are kept
#' @param v character vector of nodes over which to calculate network potential.
#' @param neighbors named list containing the neighbors for each node of graph g. If not provided,
#'        it will be computed
#' @return dataframe containing network potential for each of the inputed gene names.
#'

calc_np_all_legacy <- function(exp, g, v = as.character(names(igraph::V(g))), neighbors = NULL) {

  #first add expression to the subgraph.
  g <- add_expression(exp = exp, g = g)

  vertices <- as.character(names(igraph::V(g))) #in most cases this will be the same as `v`

  # remove any names of exp that are not in the graph.
  exp <- exp[names(exp) %in% vertices]

  # need to do the same thing for v (because sometimes there will be neighbors of a given node that aren't in the ppi)
  v <- v[v %in% vertices]

  #get a list of neighbors for each node
  if(is.null(neighbors)) {
    neighbors <-
      lapply(v,
             get_neighbors, g = g)
    names(neighbors) <- v

  }
  #loop over all vertices
  np_vec <- vector(mode = "numeric", length = length(v))
  names(np_vec) <- v
  vertex_list <- igraph::V(g) #slightly different than "vertices" - used for indexing
  for(i in v) {
    neighbors_named <- as.character(names(vertex_list[neighbors[[i]]])) #grab named vec of neighbors for each vertex
    c_j <- sum(exp[neighbors_named]) #sum up the concentration of all neighbors
    c_i <- exp[i]
    np_vec[i] <- calc_np(c_i = c_i, c_j = c_j)
  }

  #gonna do some data wrangling to get the vertices in the same order as the input expression vector
  if(length(np_vec) < length(exp)) {
    exp <- exp[names(exp) %in% names(np_vec)]
  }
  np_vec <-np_vec[names(exp)]

  return(np_vec)
}

#' function to calculate the network potential for each protein in a user-provided vector - cpp internal version
#'
#' @param exp expression vector - assumed to be a named vector where the values are expression and the names are the gene name
#' @param g igraph object - will be filtered so that only nodes found in both exp and g are kept
#' @param v character vector of nodes over which to calculate network potential.
#' @param neighbors named list containing the neighbors for each node of graph g. If not provided,
#'        it will be computed
#' @return dataframe containing network potential for each of the inputed gene names.
#'
#' @export

calc_np_all <- function(exp, g, v ="default",
                         neighbors = NULL) {
  #define list of vertices to calculate np for.
  if (length(v) <= 1) {
    if(v == "default") {
      v <- as.character(names(igraph::V(g)))
      if(length(v) < 1) {
        v <- as.character(1:length(igraph::V(g)))
      }
    }
  }


  if(is.null(names(exp))) {
    names(exp) <- as.character(1:length(exp))
  }

  #Define the total list of vertices
  vertices <- as.character(names(igraph::V(g)))
  if(length(vertices) < 1) {
    vertices <- as.character(1:length(igraph::V(g)))
  }

  #first add expression to the subgraph.
  g <- add_expression(exp = exp, g = g)

  # remove any names of exp that are not in the graph.
  exp <- exp[names(exp) %in% vertices]

  # need to do the same thing for v (because sometimes there will be neighbors of a given node that aren't in the ppi)
  # actually how is that even possible? leaving it in for now
  v <- v[v %in% vertices]
  #v we will just use as a logical flag

  #re-order exp to have the same order as vertices
  exp <- exp[vertices]

  if(is.null(neighbors)) {
    neighbors <-
      lapply(vertices,
             get_neighbors, g = g)
    names(neighbors) <- vertices
  } else {
    neighbors <- neighbors[vertices] #index to make sure we don't have any list items we don't need and to make sure they are in the same order as exp and vertices
  }

  #run cpp function to do the actual calculation on each node
  np_vec <- fcalc_np_all(neighbors = neighbors, vertices = vertices, v = v, exp = exp)
  #gonna do some data wrangling to get the vertices in the same order as the input expression vector
  # if(length(np_vec) < length(exp)) {
  #   exp <- exp[names(exp) %in% names(np_vec)]
  # }

  return(np_vec)
}

#' calculate network potential for one node.
#'
#' @param c_i expression for a given node.
#' @param c_j vector of expressions for each neighbor of c_i
#'
#' @export

calc_np <- function(c_i, c_j) {
  g_i <- c_i*log(c_i/(sum(c_j)))
  return(g_i)
}

