#' Convert from most other representations of gene name to gene.symbol
#'
#' @param x vector of ensemble.gene ids, ensemble.peptide ids, ensemble.transcript
#'      ids or entrez gene ids
#' @param edb ensemble database object
#'
#' @return vector of gene symbols
#'
#' @examples
#'\donttest{
#' #1) from numeric formatted entrez id
#' as_gene_symbol(1956)
#' #2) from character formatted entrez id
#' as_gene_symbol("1956")
#' #3) from ensemble gene id
#' as_gene_symbol("ENSG00000146648")
#' #4) From a vector of entrez ids
#' as_gene_symbol(c("123", "1956", "2012"))
#'}
#' @export

as_gene_symbol <- function(x, edb = NULL) {

#users can pass an edb object. If they don't - it uses the one below
  if(is.null(edb)) {
    edb <- EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79
  }
  #in case the input is a list - this will coerce to a vector
  if(is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  key_type = detect_inputtype(x)
  #Sometimes ensemble ids have a leading numeric string denoting some version number..
  #lets remove that
  if(is_ensembl(x)){
    x <- stringr::str_replace(x, "^[:digit:]*\\.", "")
  }

  #If its an entrez_id it may be numeric in which case we should coerce to character
  if(is.numeric(x)){
    x <- as.character(x)
  }
  #connect to ensemble db to retrieve annotations
  geneIDs <- ensembldb::select(edb,
                               keys = x,
                               keytype = key_type,
                               columns = c(key_type,"SYMBOL"))

  #store user-supplied vector as a tibble to join. - this step just makes sure that the keys are in the correct order.
  x_tibble <- tibble::tibble(x)
  colnames(x_tibble) <- key_type

  #if entrez- need to coerce back to a number for joining
  if(key_type == "ENTREZID") {
    x_tibble$ENTREZID <- as.numeric(x_tibble$ENTREZID)
  }


  gene_ids<- dplyr::left_join(x_tibble, geneIDs)
  y = gene_ids$SYMBOL #grab vector
  return(y)
}

#' Determine which format of gene is used to specify by user-defined seed proteins
#'
#' @param x vector of gene symbols
#'
#' @return "gene_symbol", "entrez_id", "ensemble_id" or "other"

detect_inputtype <- function(x) {
  if(is_ensembl(x)) { #still need to sort out whether its proteins or genes.
    return(ensembl_type(x))
  } else if (is_entrez(x)) {
    return("ENTREZID")
  }
}

#' Determine if a character vector contains ensembl gene_ids
#'
#' @param x vector or single gene symbol
#'
#' @return logical
is_ensembl <- function(x) {
  if(!is.character(x)){
    return(FALSE)
  }
  if(all(stringr::str_detect(x, "ENS.0000"))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Determine if ensembl id is a Protein, gene, or transcript_id
#'
#' @param x vector or single gene symbol
#'
#' @return character: "PROTEINID", "GENEID", "TRANSCRIPTID"

ensembl_type <- function(x) {
  if(any(stringr::str_detect(x, "(?<=ENS)P"))) {
    return("PROTEINID")
  } else if(any(stringr::str_detect(x, "(?<=ENS)G"))) {
    return("GENEID")
  } else if(any(stringr::str_detect(x, "(?<=ENS)T"))) {
    return("TRANSCRIPTID")
  } else {
    stop("input must be an ensemble gene, peptide, or transcript id")
  }
}

#' Determine if a character vector contains entrez gene_ids
#'
#' @param x vector or single gene symbol
#'
#' @return logical

is_entrez <- function(x) {

  if(is.character(x)){
    return(any(!stringr::str_detect(x, "[:alpha:]"))) #if its a character vector with no alphabetical characters - assume the user wants them to be entrez genes.
  }
  if(is.numeric(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}




