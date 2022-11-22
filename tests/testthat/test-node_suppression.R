load(system.file("test_data/toy_graph.Rda", package = "disruptr"))
#np_all
seeds <- c("OLR1", "APP", "VAV2", "ITGAV", "JAG1", "APOH")
toy_exp <- c(4.9, 9.9, 1.0, 0.2, 7.5, 8.4)
names(toy_exp) <- seeds

expected_np_old <- c(15.7, -4.8, Inf, -0.94, -2.23, -1.38)
names(expected_np_old) <- seeds
v_check <- seeds[5]
g <- igraph::induced_subgraph(g, seeds)

np_diff <- node_repression(g, v_rm = v_check, exp = toy_exp,
                           state_function = calc_np_all)
np_diff <- np_diff[,1]
expected_np_diff <- c(0.0, 6.21, 0.0, 0.082, 2.23, 0.0)
names(expected_np_diff) <- seeds

expected_np_diff <- expected_np_diff[names(np_diff)]


test_that("node_repression produces expected results", {
  expect_equal(np_diff,
               expected_np_diff, tolerance = 0.1)

})

######Multiple Nodes#####

np_diff <- node_repression(g, v_rm = seeds[1:3], exp = toy_exp,
                           state_function = calc_np_all)

test_that("removing multiple nodes works", {
  expect_equal(ncol(np_diff), 3)
  expect_equal(nrow(np_diff), 6)
})

####Every Node#####

np_diff <- node_repression(g, v_rm = seeds, exp = toy_exp, state_function = calc_np_all)

test_that("removing every node works", {
  expect_equal(ncol(np_diff), 6)
  expect_equal(nrow(np_diff), 6)
})

