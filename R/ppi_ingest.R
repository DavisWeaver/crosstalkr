#' Prepare Stringdb for use in analyses
#'
#' Basically a wrapper around the get_graph method from the stringdb package
#'
#' @param cache A filepath to a folder downloaded files should be stored
#' @param edb ensemble database object
#' @param min_score minimum connectivity score for each edge in the network.
#' @param version stringdb version
#' @param species species code either using latin species name or taxon id
#' @return igraph object built from the adjacency matrix downloaded from stringdb.
#'
#' @importFrom rlang .data
#' @export

prep_stringdb <- function(cache = NULL,
                          edb = "default",
                          min_score = 0,
                          version = "11.5", species = "homo sapiens"){
  #clean up params
  if(is.numeric(version)) {version <- as.character(version)}
  #if they provide a character version of taxon id just convert to numeric
  if(is.character(species) & !grepl("\\D", species)) {species <- as.numeric(species)}
  #if they didn't provide an edb object
  if(edb == "default") {
    edb <- EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79
  }
  if(!is.numeric(species)) {
    species = to_taxon_id(species)
  }

  if(!file.exists(paste0(cache, "/", species, "stringdb.Rda")) &
     !file.exists(paste0(cache, species, "stringdb.Rda"))) {

    message(paste0("Downloading stringdb ", species, " v", version))

    df <- STRINGdb::STRINGdb$new(version = version, species = species,
                                 score_threshold = min_score)
    g <- try(df$get_graph())
    if(inherits(df, "try-error")) {
      stop("unable to download stringdb, please try again later")
    }


    #If messing about with humans convert ensemble_ids to gene_ids and then re-converting to graph
    if(species == "homo sapiens" | species == 9606) {
      df <- as.data.frame(igraph::get.edgelist(g))
      colnames(df) <- c("protein1", "protein2")
      message("converting ensemble_ids to gene_ids")
      #Lets convert the ensemble_ids to gene_ids.
      df <- dplyr::mutate(df,
                          protein1 = as_gene_symbol(x = .data$protein1, edb = edb),
                          protein2 = as_gene_symbol(x = .data$protein2, edb = edb))
      #depending on the EDB used you end up with some NAs that you need to clear out
      df <- dplyr::filter(df, !is.na(.data$protein1), !is.na(.data$protein2))

      message("converting to igraph object")
      #Convert to igraph object
      g  <- igraph::graph_from_data_frame(df, directed = FALSE)
      g <-
        igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

    }


    if(!is.null(cache)) {
      save(g, file = paste0(cache, "/", species, "stringdb.Rda"))
    }
  } else {
    message("using cached version of stringdb Homo Sapeins v11.0")
    load(file = paste0(cache, "/", species, "stringdb.Rda"))
  }
  return(g)
}

#' Prepare biogrid for use in analyses
#'
#' @inheritParams prep_stringdb
#'
#' @return igraph object built from the adjacency matrix downloaded from thebiogrid.org.
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
#' Only works with the human PPI. min_score parameter only applies to strindb
#'
#' @inheritParams prep_stringdb
#'
#' @return igraph object corresponding to PPI following union
#' @export

ppi_union <- function(cache = NULL, min_score = 0, edb = "default") {
  g_biogrid <- prep_biogrid(cache = cache)
  g_string <- prep_stringdb(cache = cache, edb = edb, min_score = min_score)
  g <- igraph::union(g_biogrid, g_string)
  return(g)
}

#' Function to allow users to choose the intersection of stringdb and biogrid
#' Only works with the human PPI. min_score parameter only applies to strindb
#'
#' @inheritParams prep_stringdb
#'
#' @return igraph object corresponding to PPI following intersection
#'
#' @export

ppi_intersection <- function(cache = NULL, min_score = 0, edb = "default") {
  g_biogrid <- prep_biogrid(cache = cache)
  g_string <- prep_stringdb(cache = cache, edb = edb, min_score = min_score)
  g <- igraph::intersection(g_biogrid, g_string)
  return(g)
}

#' Helper function to load requested PPI w/ parameters
#'
#' @inheritParams prep_stringdb
#'
#' @param union bool
#' @param intersection bool
#' @param ppi str
#'

load_ppi <- function(cache, union = FALSE, intersection = FALSE, species = "9606", min_score, ppi= "stringdb") {
  if(union & (tolower(species) == "homo sapiens" | as.character(species) == "9606")) {
    g <- ppi_union(cache = cache, min_score = min_score)
  } else if(intersection & (tolower(species) == "homo sapiens" | as.character(species) == "9606")) {
    g <- ppi_intersection(cache = cache, min_score = min_score)
  } else if(ppi == "biogrid" & (tolower(species) == "homo sapiens" | as.character(species) == "9606")) { #first 3 options are only feasible if the species is human
    g <- prep_biogrid(cache = cache)
  } else if (ppi == "stringdb") {
    g <- prep_stringdb(cache = cache, min_score = min_score, species = species)
  } else {
    stop("ppi must be either 'biogrid' or 'stringdb'")
  }
  return(g)
}



#' helper to convert user-inputs to ncbi reference taxonomy.
#'
#' @param species user-inputted species
#'
#' @importFrom magrittr %>%
#' @return string corresponding to taxon id
#' @export

to_taxon_id <- function(species) {

  species <- stringr::str_to_lower(species) #keep it consistent

  #download reference data from string
  reference_df <- supported_species()

  #select the taxon id for that species.
  taxon_id <- dplyr::filter(reference_df,
                            .data$string_name == species |
                              .data$ncbi_name == species) %>%
    dplyr::rename(taxon_id = "#taxon_id") %>%
    dplyr::select(.data$taxon_id) %>%
    unlist()

  if(is.na(taxon_id[1])) {
    stop("Invalid species specification, for a list of supported species, call crosstalkr::supported_species()")
  }
  return(unname(taxon_id))

}


#' returns a dataframe with information on supported species
#'
#' @return dataframe
#' @importFrom magrittr %>%
#'
#' @export

supported_species <- function() {

  if(file.exists(system.file("species_reference.Rda", package = "crosstalkr"))) {
    load(system.file("species_reference.Rda", package = "crosstalkr"))
  } else {
    df <- readr::read_tsv(file = "https://stringdb-static.org/download/species.v11.5.txt")
    df <- dplyr::filter(df, .data$STRING_type == "core")
    reference_df <- dplyr::rename(df, string_name = .data$STRING_name_compact,
                                  ncbi_name = .data$official_name_NCBI) %>%
      dplyr::mutate(string_name = stringr::str_to_lower(.data$string_name),
                    ncbi_name = stringr::str_to_lower(.data$ncbi_name))
  }

  return(reference_df)

}





