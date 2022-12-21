#Lets take a look at a smaller graph
v1 = c(1,1,1,0,0,0,0,0)
v2 = c(0,0,0,1,1,0,0,0)
v3 = c(1,1,1,0,1,0,0,1)
v4 = c(0,0,0,1,1,1,1,1)
v5 = c(0,0,0,1,0,0,0,1)
v6 = c(0,1,0,1,0,1,0,1)
v7 = c(1,0,1,0,1,0,1,0)
v8 = c(1,1,1,1,0,0,1,1)
w = matrix(data = c(v1,v2,v3,v4,v5,v6,v7,v8), ncol = 8, nrow = 8)

g <- igraph::graph_from_adjacency_matrix(w)
seeds <- c(1,3,5)

test_that("Matched seeds are the same length as input seeds", {
  expect_equal(length(match_seeds(g = g, seed_proteins = seeds, n = 1)[[1]]), 3)
})

#setup test to compare distribution of seeds to what we expect
t1 = unlist(match_seeds(g = g, seed_proteins = c(1), n = 10000))

#seeds 2 3 5 6 and 7 should have been possible selections - each with a uniform chance of being selected, all others 0
test_that("Matched seeds are within 2 degree of a given seed", {
  for(i in 1:igraph::vcount(g)) {
    if(i %in% c(1,4,8)) {
      expect_equal(sum(t1 == i)/ length(t1), 0, tolerance = 0)
    } else {
      expect_equal(sum(t1== i)/ length(t1), 0.2, tolerance = 0.02)
    }
  }
})

#lets take a look at what happens on larger graphs - this is a weird case where the degree is pretty
#much always the same- turns out to be kind of challenging.

g <- igraph::sample_gnp(n = 1000, p = 100/1000)

test_that("Matched seeds are the same length as input seeds", {
  expect_equal(length(match_seeds(g = g, seed_proteins = c(1,3,5), n = 1)[[1]]), 3)
  expect_equal(length(match_seeds(g = g, seed_proteins = c(12,14), n = 1)[[1]]), 2)
  expect_equal(length(match_seeds(g = g, seed_proteins = 50, n = 1)[[1]]), 1)
  expect_equal(length(match_seeds(g = g, seed_proteins = c(21,22,23,24,25), n = 1)[[1]]), 5)
})

#Lets test the parent function
#bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 10)
test_that("bootstrap_null runs without errors",  {
  expect_true(is.data.frame(bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 10)[[1]]))
  expect_message(bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 10))
  expect_true(is.data.frame(bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 100, agg_int = 10)[[1]]))
})

test_that("bootstrap_null runs when you provide invalid vertex ids", {
  expect_true(is.data.frame(bootstrap_null(seed_proteins = c(1,3,5,9,12,15,1002), g = g, n = 10)[[1]]))
})

#Lets see what happens when you use an actual ppi (albeit a tiny one)

# load(system.file("test_data/toy_graph.Rda", package = "crosstalkr"))
# vertices <- names(igraph::V(g))
# seeds <- vertices[c(1,6,12)]
#
# test_seeds <- unlist(match_seeds(g=g, seed_proteins = seeds, n =1000))
#
# test_that("matched seeds have zero probability of being the original seeds", {
#   for(i in seeds) {
#     expect_equal(sum(test_seeds == i)/ length(test_seeds), 0, tolerance = 0)
#   }
# })
