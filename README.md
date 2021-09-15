# crosstalkr

<!-- badges: start -->
[![R-CMD-check](https://github.com/DavisWeaver/crosstalkr/workflows/R-CMD-check/badge.svg)](https://github.com/DavisWeaver/crosstalkr/actions)
[![Codecov test coverage](https://codecov.io/gh/DavisWeaver/crosstalkr/branch/main/graph/badge.svg)](https://codecov.io/gh/DavisWeaver/crosstalkr?branch=main)
<!-- badges: end -->

R package for the identification of functionally important subnetworks 

Crosstalkr provides a general implementation of a random-walk with restarts on graph structured data. 
We also provide user-friendly implementations of the common use-case of using random-walk with restarts to identify subnetworks of biological protein-protein interaction databases. 
Given a user-defined set of seed proteins, the main `compute_crosstalk` function will compute affinity scores for all other proteins in the network. 
It will then compute a null distribution using a permutation test and compare the computed affinity scores to the null distribution to identify proteins with a statistically signficant association to the user-defined seed-proteins.
This user-friendly implementation currently only supports human protein-protein interaction networks. 
Soon, users will be able to select from a list of additional supported species. 


# Use

To install, use the following code: 

```
install.packages("remotes") #can skip if already installed 
remotes::install_github("https://github.com/DavisWeaver/crosstalkr")
```

Given a set of user-provided set of seed proteins, crosstalkr will identify enriched subnetworks of genes that have a high affinity for the provided seeds. 

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


More detailed methodology can be found here: "Nibbe et al., An Integrative -omics Approach to Identify Functional Sub-Networks in Human Colorectal Cancer, Comput Biol 6(1): e1000639. doi:10.1371/journal.pcbi.1000639"

# Contact

Please feel free to submit issues here. You can also contact me at davis.weaver@case.edu if you have any questions. 


