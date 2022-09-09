---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: '"Crosstalkr: an R package for the identification of related nodes in biological networks"'
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
date: "September, 07"
year: "2022"
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

Crosstalkr is eesigned to facilitate the identification of functional subnetworks in human or non-human interactomes. 
Given a set of user-provided set of seed proteins, crosstalkr will identify enriched subnetworks of proteins that have a high affinity for the provided seeds. 
It is a free, open-source R package designed to allow users to integrate functional analysis using the protein-protein interaction network into existing bioinformatic pipelines. 
This is accomplished using random walks with restarts, starting at the user-provided seed proteins. 
Random walks are implemented using sparse matrix multiplication to facilitate fast execution. 
Affinity scores from a given random walk with restarts are compared to a bootstrapped null distribution to assess statistical significance. 
The default behavior evaluates the human interactome to identify functionally important subnetworks given a set of user-defined seed proteins. 
However, users can also provide a different graph, allowing for flexible evaluation of graph or network-structured data. 
Further, users can evaluate more than 1000 non-human protein-protein interaction networks thanks to integration with StringDB.
Crosstalkr is an extension of the crosstalker webapp developed by Neoproteomics and Case Western Reserve University School of Medicine. 

# Statement of Need

# Design and Data Sources

# Functionality

# Outlook

# Acknowledgements

We acknowledge contributions from Mark Chance and Mehmet Koyuturk.

# References

# Citations

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"
