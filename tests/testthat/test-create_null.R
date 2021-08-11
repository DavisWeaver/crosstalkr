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

test_seeds <- unlist(match_seeds(g=g, seed_proteins = seeds, n =1000))
test_that("Matched seeds are randomly selected if igraph::vcount(g) < 50", {
  expect_equal(sum(test_seeds == 1)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 2)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 3)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 4)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 5)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 6)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 7)/ length(test_seeds), 1/8, tolerance = 0.05)
  expect_equal(sum(test_seeds == 8)/ length(test_seeds), 1/8, tolerance = 0.05)
})

#lets take a look at what happens on larger graphs - this is a weird case where the degree is pretty
#much always the same- turns out to be kind of challenging.

g <- igraph::sample_gnp(n = 1000, p = 10/1000)

test_that("Matched seeds are the same length as input seeds", {
  expect_equal(length(match_seeds(g = g, seed_proteins = c(1,3,5), n = 1)[[1]]), 3)
  expect_equal(length(match_seeds(g = g, seed_proteins = c(12,14), n = 1)[[1]]), 2)
  expect_equal(length(match_seeds(g = g, seed_proteins = 50, n = 1)[[1]]), 1) #this one breaks still
  expect_equal(length(match_seeds(g = g, seed_proteins = c(21,22,23,24,25), n = 1)[[1]]), 5)
})

seeds <- 20:40
ds <- igraph::degree_distribution(g, v = seeds)
t1 <- match_seeds(g, seed_proteins = seeds, n = 100)
# this code will compute a degree distribution from the computed seeds
t2 <- list()
for(i in 1:length(t1)) {
  t2[[i]] <- igraph::degree(g=g,v = t1[[i]])
}
seeds_props <- table(unlist(t2))/ length(unlist(t2))


test_that("matched seeds actually track the input distribution", {
  expect_equal(min(igraph::degree(g, seeds)), as.numeric(names(seeds_props[1])))
  expect_equal(max(igraph::degree(g, seeds)), as.numeric(names(seeds_props[length(seeds_props)])))

})


#Lets test the parent function
bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 10)
test_that("bootstrap_null runs without errors",  {
  expect_true(is.data.frame(bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 10)[[1]]))
  expect_message(bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 10))
  expect_true(is.data.frame(bootstrap_null(seed_proteins = c(1,3,5,9,12,15), g = g, n = 100, agg_int = 10)[[1]]))
})




