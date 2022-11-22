load(system.file("test_data/toy_graph.Rda", package = "disruptr"))

#np_all
seeds <- c("OLR1", "APP", "VAV2", "ITGAV", "JAG1", "APOH")
toy_exp <- c(4.9, 9.9, 1.0, 0.2, 7.5, 8.4)
names(toy_exp) <- seeds

expected_np <- c(15.7, -4.8, Inf, -0.94, -2.23, -1.38)
names(expected_np) <- seeds

g <- igraph::induced_subgraph(g, seeds)

testnp <- calc_np_all(g = g, exp = toy_exp)
testnp <- testnp[names(expected_np)]

test_that("calc_np_all returns expected output", {
 expect_lte(length(calc_np_all(toy_exp, g)), length(toy_exp))
 expect_true(is.numeric(calc_np_all(toy_exp,g)))
 expect_equal(testnp, expected_np, tolerance = 0.1)
})


g_new <- igraph::delete.vertices(g, "JAG1")
seeds_new <- seeds[seeds != "JAG1"]
expected_np_new <- c(15.67, 1.39, Inf, -0.86, -1.38)
names(expected_np_new) <- seeds_new
testnp <- calc_np_all(g = g_new, exp = toy_exp)
testnp <- testnp[names(expected_np_new)]

test_that("calc_np_all works with different length arguments", {
expect_equal(testnp, expected_np_new, tolerance = 0.1)
})

expected_np <- c(15.67, -4.81, Inf,0,0,0)
names(expected_np) <- seeds
testnp <- calc_np_all(g = g, exp = toy_exp, v = c("OLR1", "APP", "VAV2"))
testnp <- testnp[names(expected_np)]
test_that("calc_np_all can process subsets of a graph", {
  expect_equal(testnp,
               expected = expected_np, tolerance = 0.1)
})

test_that("calc_np returns expected output", {
  expect_equal(calc_np(c_i= 10, c_j = 20), -6.931472, tolerance = 0.0001)
  expect_gte(calc_np(c_i = 10, c_j = 5), 0)
  expect_lte(calc_np(c_i = 10, c_j = 11), 0)
  expect_true(is.infinite(calc_np(c_i = 1, c_j = 0)))
})




