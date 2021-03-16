#bring in functions
load_all()
library(magrittr)
library(dplyr)

g <- prep_stringdb(cache = "test_data")
x <- unlist(igraph::vertex.attributes(g), use.names = FALSE)

#Lets use one of the hallmarks datasets or something
h_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")

seed_proteins <- h_gene_sets %>% filter(gs_name == "HALLMARK_ANGIOGENESIS") %>%
  select(gene_symbol) %>% unlist()
save(seed_proteins, file = "./test_data/seed_proteins.Rda")


##Also want to try random seeds

random_seeds <- sample(x, size = 100)
save(random_seeds, file = "./test_data/random_seeds.Rda")

