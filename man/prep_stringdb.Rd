% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppi_ingest.R
\name{prep_stringdb}
\alias{prep_stringdb}
\title{Prepare Stringdb for use in analyses}
\usage{
prep_stringdb(
  cache = NULL,
  edb = "default",
  min_score = 200,
  version = "11.5",
  species = "homo sapiens"
)
}
\arguments{
\item{cache}{A filepath to a folder downloaded files should be stored}

\item{edb}{ensemble database object}

\item{min_score}{minimum connectivity score for each edge in the network.}

\item{version}{stringdb version}

\item{species}{species code either using latin species name or taxon id}
}
\value{
igraph object built from the adjacency matrix downloaded from stringdb.
}
\description{
Basically a wrapper around the get_graph method from the stringdb package
}
