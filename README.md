# crosstalkr
R package for the identification of functionally important subnetworks 

<!-- badges: start -->
[![R-CMD-check](https://github.com/DavisWeaver/crosstalkr/workflows/R-CMD-check/badge.svg)](https://github.com/DavisWeaver/crosstalkr/actions)
<!-- badges: end -->

Crosstalkr is an extension of the crosstalker webapp developed by Neoproteomics and Case Western Reserve University School of Medicine. 
It is a free, open-source R package designed to provide the same functionality as the crosstalker webapp in a modular fashion that can be incorporated into existing bioinformatic pipelines. 

# Use

Given a set of user-provided set of seed proteins, crosstalkr will identify enriched subnetworks of genes that have a high affinity for the provided seeds. 
This is accomplished using random walks with restarts, starting at the user-provided seed proteins. 
Random walks are implemented using sparse matrix multiplication to facilitate fast execution. 

Affinity scores from a given random walk with restarts are compared to a bootstrapped null distribution to assess statistical significance. 
More detailed methodology can be found here: "Nibbe et al., An Integrative -omics Approach to Identify Functional Sub-Networks in Human Colorectal Cancer, Comput Biol 6(1): e1000639. doi:10.1371/journal.pcbi.1000639"

# Contact

Please feel free to submit issues here. You can also contact me at davis.weaver@case.edu if you have any questions. 


