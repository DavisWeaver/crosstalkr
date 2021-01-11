#' Convert from ensembl.gene to gene.symbol
#'
#' @param x vector of ensemble.gene ids
#' @param edb ensemble database object
#'
#' @return vector of gene.symbol

ensembl_convert <- function(x, edb = EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79) {

  if(is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }

  #Sometimes ensemble ids have a leading numeric string denoting some version number..
  #lets remove that

  x <- stringr::str_replace(x, "^[:digit:]*\\.", "")

  #connect to ensemble db to retrieve annotations
  geneIDs <- ensembldb::select(edb,
                               keys = x,
                               keytype = "PROTEINID",
                               columns = c("PROTEINID","SYMBOL","GENEID"))
  gene_ids<- dplyr::left_join(tibble::tibble(PROTEINID = x), geneIDs) #grab vector
  y = gene_ids$SYMBOL
  return(y)
}


