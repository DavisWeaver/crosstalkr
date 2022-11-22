#' Function to eliminate a node from a network g and calculate the change in some measure of network state
#'
#' this function is still under development.
#'
#' @param g igraph network object
#' @param v_rm index of vertices to remove
#' @param exp expression vector for nodes in graph g
#' @param state_function function to use to calculate network state before and after node_repression
#' @param neighbors_only logical designating whether state function should be calculated for all nodes or just neighbors
#' @param ... additional parameters passed to state function.
#' @export

node_repression <- function(g, v_rm, exp, state_function = calc_np_all,
                            neighbors_only = TRUE, ...) {


  #add expression now to fix indexing
  g <- add_expression(exp=exp, g=g)
  #grab all vertices of g
  vertices <- as.character(names(igraph::V(g)))

  #limit v_rm to those present in vertices
  v_rm <- v_rm[v_rm %in% vertices]

  #Get a list of neighbors in the graph - this will speed things up a lot if we are suppressing tons of nodes simultaneously
  neighbors_all <- lapply(vertices,
                          get_neighbors, g = g)
  names(neighbors_all) <- vertices

  out_mat <- Matrix::Matrix(0, nrow = length(igraph::V(g)), ncol = length(v_rm),
                            sparse = TRUE)
  row.names(out_mat) <- vertices
  colnames(out_mat) <- v_rm

  for(i in 1:length(v_rm)) {
    #define v to be either all nodes or neighbors of v_rm
    if(neighbors_only == TRUE) {
      v <- neighbors_all[[v_rm[i]]] #get neighbors
      v <- as.character(names(igraph::V(g)[v]))
    } else {
      v <- as.character(names(igraph::V(g)))
    }
    #g_new <- igraph::delete_vertices(g = g, v = v_rm[i])
    exp_new <- exp
    exp_new[v_rm[i]] <- 0
    #calculate state vector for original graph
    s_old <- state_function(exp = exp, g = g, v=c(v,v_rm[i]), neighbors = neighbors_all) #calculate network potential for the deleted note in addition to the neighbors

    #create a 0 value deleted node to add on
    # delete_node <- 0
    # names(delete_node) <- v_rm[i]
    s_new <- state_function(exp = exp_new, g = g, v = v, neighbors = neighbors_all) #calc new network potential for the neighbors + 0 for the deleted node

    s_new <- s_new[names(s_old)]
    s_diff <- s_new-s_old #calculate the diff and save
    out_mat[,i] <- s_diff

  }
  #convert all infinites to 0
  out_mat[is.infinite(out_mat)] <- 0
  return(out_mat)
}
