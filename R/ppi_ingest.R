#' Prepare Stringdb for use in analyses
#'
#'
#'
#' @param cache A filepath to a folder downloaded files should be stored, inherits from user-available functions
#' @param edb ensemble database object
#' @param min_score minimum connectivity score for each edge in the network.
#' @return list containing Adjacency matrix from stringdb dataset and igraph object built from the adjacency matrix.
#'
#' @importFrom rlang .data
#' @export

prep_stringdb <- function(cache = NULL,
                          edb = "default",
                          min_score = NULL){
  #if they didn't provide an edb object
  if(edb == "default") {
    edb <- EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79
  }
  if(!file.exists(paste0(cache, "/stringdb.Rda"))) {
    message("Downloading stringdb Homo Sapiens v11.0")
    df <-try(readr::read_delim("https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz",
                            delim = " "))

    if(inherits(df, "try-error")) {
      stop("unable to download stringdb, please try again later")
    }

    message("converting ensemble_ids to gene_ids")
    #Lets convert the ensemble_ids to gene_ids.
    df <- dplyr::mutate(df,
                        protein1 = as_gene_symbol(x = .data$protein1, edb = edb),
                        protein2 = as_gene_symbol(x = .data$protein2, edb = edb))

    message("filtering proteins with a certain min_score")
    #filter out nodes below a given min score
    if(is.numeric(min_score)) {
      df <- dplyr::filter(df, .data$combined_score > min_score)
    }

    message("converting to igraph object")
    #Convert to igraph object
    g  <- igraph::graph_from_data_frame(df, directed = FALSE)
    g <-
      igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

    if(!is.null(cache)) {
      save(g, file = paste0(cache, "/stringdb.Rda"))
    }
  } else {
    message("using cached version of stringdb Homo Sapeins v11.0")
    load(file = paste0(cache, "/stringdb.Rda"))
  }
  return(g)
}

#' Prepare biogrid for use in analyses
#'
#' @inheritParams prep_stringdb
#'
#' @return list containing Adjacency matrix from stringdb dataset and igraph object built from the adjacency matrix.
#'
#' @export

prep_biogrid <- function(cache = NULL) {

  if(!file.exists(paste0(cache, "/biogrid.Rda"))) {
    tmp <- tempdir()

    #Download most recent version of the biogrid
    message("Downloading biogrid version 3.5.171")
    download.file("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.171/BIOGRID-ORGANISM-3.5.171.tab2.zip",
                  destfile = paste0(tmp, "/biogrid.zip"))

    #Unzip only the homosapiens portion of the biogrid zip file and delete big zip file
    unzip(paste0(tmp, '/biogrid.zip'),
          files = c("BIOGRID-ORGANISM-Homo_sapiens-3.5.171.tab2.txt"),
          exdir = paste0(tmp, "/unzip"))
    file.remove(paste0(tmp, "/biogrid.zip"))

    #Read in biogrid file from temp directory.
    biogrid <-
      try(read.delim(
        paste0(tmp, "/unzip/BIOGRID-ORGANISM-Homo_sapiens-3.5.171.tab2.txt"),
        header = TRUE
      ))
    if(inherits(biogrid, "try-error")) {
      stop("error downloading the biogrid protein-protein interaction database. Please try again later. ")
    }
    biogrid <-
      biogrid[, 8:9] # isolate Official Symbol Interactor A & B columns

    #convert to graph class and simplify
    g  <- igraph::graph_from_data_frame(biogrid, directed = FALSE)
    g <-
      igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

    if(!is.null(cache)) {
      save(g, file = paste0(cache, "/biogrid.Rda"))
    }

  } else {
    message("using cached version of biogrid v3.5.171")
    load(file = paste0(cache, "/biogrid.Rda"))
  }
  return(g)
}

#' Function to allow users to choose the union of stringdb and biogrid
#'
#'

#' Helper function for first-time use of crosstalkr package
#'
#' @inheritParams prep_stringdb
#'
#' @return directory on users computer containing the different adjacency matrices for future use.
#'
#' @export
#'
setup_init <- function(cache = NULL, min_score) {

  #Functons are written to return a tibble - this use will ensure a df is not printed
  tmp_var1 <- prep_biogrid(cache = cache)
  tmp_var2 <- prep_stringdb(cache = cache, min_score = min_score)
}


