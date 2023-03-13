# crosstalkr

<!-- badges: start -->
[![R-CMD-check](https://github.com/DavisWeaver/crosstalkr/workflows/R-CMD-check/badge.svg)](https://github.com/DavisWeaver/crosstalkr/actions)
[![Codecov test coverage](https://codecov.io/gh/DavisWeaver/crosstalkr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/DavisWeaver/crosstalkr?branch=main)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![CRAN version](http://www.r-pkg.org/badges/version/crosstalkr)](https://CRAN.R-project.org/package=crosstalkr)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/crosstalkr)](https://CRAN.R-project.org/package=crosstalkr)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FDavisWeaver%2Fcrosstalkr&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
<!-- badges: end -->

# R package for the analysis of biological networks. 

Crosstalkr provides a unified toolkit for drug target and disease subnetwork identification. Crosstalkr enables users to download and leverage high-quality protein-protein interaction networks from online repositories. 
Users can then filter these large networks into manageable subnetworks. 
Finally, users can perform in-silico repression experiments to assess the relative importance of nodes in their network.

## PPI ingestion and customization

Crosstalkr allows direct access to the STRINGDB and Biogrid PPI resources. Thanks to integration with stringdb, users can evaluate biological networks from 1540 different species. For a list of supported species - just call:

```
crosstalkr::supported_species()
```

Users can also take the intersection or union of the stringdb and biogrid PPIs. 

## Graph Filtering Methods

Crosstalkr faciliates graph reduction based on node value ranks. 
Node values can be provided by the user (as in `gfilter.value`).
Users can also specify any method found in the `igraph` package that generates node values (i.e. `igraph::degree` or `igraph::betweenness`).
We also provide a custom method for node ranking that we developed in our lab, termed network potential (`gfilter.np`)
Crosstalkr provides a general implementation of a random-walk with restarts on graph structured data. 
We also provide user-friendly implementations of the common use-case of using random-walk with restarts to identify subnetworks of biological protein-protein interaction databases. 
Given a user-defined set of seed proteins, the main `compute_crosstalk` function will compute affinity scores for all other proteins in the network. 
It will then compute a null distribution using a permutation test and compare the computed affinity scores to the null distribution to identify proteins with a statistically significant association to the user-defined seed-proteins.

## Node ranking via in-silico repression

In silico repression is implemented by the `node_repression` function. 
Users must specify a state function that scores nodes. Each node in `v_rm`, will be systematically removed from the network. 
The provide state function will be applied to re-calculate network state and then the difference in total state value (sum of all nodes) will be computed. 

See https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1008755 for more details on in-silico repression. 

## Other functionality


# Use

To install, use the following code: 

```
install.packages("crosstalkr")

```

For the latest development version:
```
install.packages("remotes") #can skip if already installed 
remotes::install_github("https://github.com/DavisWeaver/crosstalkr")
```

Given a set of user-provided set of seeds, crosstalkr will identify enriched an enriched subgraph where all nodes have a high affinity for the provided seeds. 

crosstalkr is optimized for use with the human cell signaling network. For example, running the code below will return a dataframe containing the user-provided seeds as well as all other proteins in the human protein-protein interaction network with a statistically significant association to these genes.

```
compute_crosstalk(c("EGFR", "KRAS"))
```

Users can use any other kind of graph-structured data, provided they are stored in an igraph object. For example:


```
g <- igraph::sample_gnp(n = 1000, p = 10/1000)
compute_crosstalk(c(1,3,5,8,10), g = g, use_ppi = FALSE)
```

We also provide rudimentary plotting functions to allow users to quickly assess the identified subnetworks: 

```
ct_df <- compute_crosstalk(c("EGFR", "KRAS"))
plot_ct(ct_df)
```

A more detailed overview of the available functionality is provided in the introductory vignette (under development). 

```
vignette(package = "crosstalkr")
```

Please use the provided biorxiv pre-print to cite. https://www.biorxiv.org/content/10.1101/2023.03.07.531526v1

# Contact

Please feel free to submit issues here. You can also contact me at davis.weaver@case.edu if you have any questions. 


