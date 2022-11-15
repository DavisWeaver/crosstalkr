---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'crosstalkr: an R package for the identification of related nodes in biological networks'
tags:
  - R
  - Graph Theory
  - Interactomics
authors:
  - name: Davis T. Weaver
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Jacob Scott
    affiliation: "1,2"
affiliations:
 - name: Case Western Reserve University School of Medicine
   index: 1
 - name: Cleveland Clinic Foundation
   index: 2
citation_author: Weaver et. al.
date: "November, 15"
year: "2022"
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

Crosstalkr is designed to facilitate the identification of functional subnetworks in graph-structured data. 
It is a free, open-source R package designed to allow users to integrate functional analysis using the protein-protein interaction network into existing bioinformatic pipelines. 
Given a set of user-provided seed proteins, crosstalkr will identify a group of proteins that have a high affinity for the provided seeds. 
This is accomplished using random walks with restarts, starting at the user-provided seed proteins. 
Affinity scores from a given random walk with restarts are compared to a bootstrapped null distribution to assess statistical significance. 
Random walks are implemented using sparse matrix multiplication to facilitate fast execution. 
The default behavior evaluates the human interactome.
However, users can also provide a different graph, allowing for flexible evaluation of graph or network-structured data. 
Further, users can evaluate more than 1000 non-human protein-protein interaction networks thanks to integration with StringDB.



# Statement of Need

In the last few decades, interest in graph-based analysis of biological networks has grown substantially. 
Protein-protein interaction networks are one of the most common biological networks, and represent the molecular relationships between every known protein and every other known protein.

Researchers have applied graph search and graph clustering algorithms to biological networks in an effort to derive disease-specific subgraphs [@nibbe_integrative_2010;@chitra_netmix2_2022;@pfeifer_gnn-subnet_2022] or identify potential drug targets [@weaver_network_2021;@martinez_drugnet_2015]. 
One of the most well-studied algorithms in this context is the random walk with restarts.
Random walks with restarts (RWR) have been used and adapted across disciplines and industries for applications ranging from internet search engines to drug target identification.  [@tong_fast_2006;@bianchini_inside_2005;@navarro_prophtools_2017]

There is a growing suite of tools available in R for analyzing graph-structured data [@details_igraph_2022; @gatto_using_2014], including a few R packages that implement RWR [@fang_dnet_2014;@valentini_ranks_2022]. 
These include the RANKS package, which provides tools for performing many graph-based node scoring algorithms.(@valentini_ranks_2022).
These tools require some understanding of graph data structures and ask the user to find, download, and manipulate the relevant biological networks into adjacency matrices or igraph objects.
Crosstalkr compromises some of the flexibility of RANKS to provide an optimized, streamlined interface to allow users to integrate interactomic analyses into their workflow. 
While users can interact directly with crosstalkr to perform RWR on any graph, the package is optimized to facilitate one-line implementation of an algorithm designed to identify functional subgraphs of protein-protein interaction networks (PPI). 

# Design and Data Sources

## compute_crosstalk

The main entrypoint for most users will be the compute_crosstalk function.
If users plan to search a supported protein protein interaction network, they are 
only required to provide a vector of seed proteins. 
In this situation, compute_crosstalk will:

  1. Download the requested PPI (or load it from the provided cache)
  2. Process the requested PPI into a sparse adjacency matrix.
  3. Perform a random walk with restart using the user provided seeds to generate affinity scores for every protein in the PPI.
  4. Perform many random walks with restarts from n random seeds with a matching degree distribution to generate a null distribution of affinity score.
  5. Compare the affinity scores to the null distribution to compute an adjusted p-value (using the method specified in p_adjust)
  6. Remove proteins with an adjusted p-value < significance_level 

Users can make use of caching to store processed PPIs and speed up future analyses substantially. 
Users can also make use of parallel computing by setting the ncores parameter > 1
A sample workflow demonstrating the ease of use is provided below. 
Here, we attempt to determine proteins that are functionally related to EGFR, KRAS, PI3K, and STAT3; proteins that are involved in growth signaling in cancer cells. 



```r
df <- compute_crosstalk(seed_proteins = c("EGFR", "KRAS","STAT3"), 
                        cache = "./data", seed_name = "joss_ex", n = 10000, 
                        significance_level = 0.99)
df %>%
  select(-c(Z,mean_p, var_p, nobs)) %>%
  slice_max(order_by = affinity_score, n = 5) %>% 
  knitr::kable(digits = 4)
```



|node    |seed | affinity_score| p_value| adj_p_value|
|:-------|:----|--------------:|-------:|-----------:|
|STAT3   |yes  |         0.2002|       0|       0e+00|
|EGFR    |yes  |         0.2001|       0|       0e+00|
|KRAS    |yes  |         0.2001|       0|       0e+00|
|C2orf72 |no   |         0.0041|       0|       1e-04|
|CCDC87  |no   |         0.0029|       0|       0e+00|

We also provide a convenience function to quickly plot the returned subgraph. Users can specify `prop_keep` to improve readability by only plotting the top x% of identified proteins, ranked by affinity score. 


```r
g <- prep_stringdb(cache = "./data")
crosstalkr::plot_ct(df, g=g, prop_keep = 0.4, label_prop = 0.2)
```

![Protein-protein interaction subnetwork for EGFR, KRAS, and STAT3. ](crosstalkr_joss_files/figure-latex/unnamed-chunk-2-1.pdf) 

## Other Features

In pursuit of a one-line interactomic analysis pipeline, we developed several convencience functions that users will likely find useful in other analyses. 
For example, the human protein-protein interaction network crosstalkr is able to detect and convert between entrez ids, uniprot names, and ensemble ids. 
Users can make use of the as_gene_symbol function to convert any common representation of gene identity into human-readable gene names.
In addition, users can make use of single-line functions that download and clean PPIs from either StringDB or Biogrid.
Crosstalkr also ships with a highly optimized implementation of the random-walk with restarts algorithm (sparseRWR), which users can apply to any graph-structured data. 

## Data sources 

Users can leverage two high quality PPIs through crosstalkr; StringDB and Biogrid [@oughtred_biogrid_2021; @szklarczyk_string_2021]. Users can run their analysis using either of these resources individually or they can take the union or intersection of these networks. While Biogrid only supports the human PPI, StringDB provides high quality PPIs for more than 1500 species [@szklarczyk_string_2021]. 
Crosstalkr provides a user-friendly interface for all of these species. 

# Acknowledgements

We acknowledge contributions from Mark Chance and Mehmet Koyuturk. 
We would also like to acknowledge the dependencies that enable crosstalkr [@wickham_dplyr_2022;@details_igraph_2022;@magrittr_magrittr_2022;@hester_withr_2022;@bates_matrix_2022;@wickham_readr_2022;@wickham_tidyr_2022;@daniel_foreach_2022;@daniel_doparallel_2022;@rainer_ensembldb_2019] and this paper [@xie__aut_knitr_2022;@wickham_ggplot2_2022]

# Citations

