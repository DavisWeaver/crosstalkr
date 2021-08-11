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

test_that("compute_crosstalk doesn't break for very small graphs", {
  expect_message(compute_crosstalk(c(1,3), g = g, use_ppi=FALSE))
  expect_true(is.data.frame(compute_crosstalk(c(1,3), g = g, use_ppi=FALSE, n=100)))
  expect_true(is.data.frame(compute_crosstalk(c(1), g = g, use_ppi=FALSE, n = 100)))
  expect_null(check_crosstalk(compute_crosstalk(c(1), g = g, use_ppi=FALSE, n = 100)))
})

g <- igraph::sample_gnp(n = 1000, p = 10/1000)

test_that("compute_crosstalk identifies crosstalkers for larger graphs", {
  expect_true(is.data.frame(compute_crosstalk(c(1,3,5,8,10), g = g, use_ppi = FALSE, n = 100)))
  expect_null(check_crosstalk(compute_crosstalk(c(1), g = g, use_ppi=FALSE, n = 100)))
  expect_true(nrow(compute_crosstalk(c(1,3,5,8,10), g = g, use_ppi = FALSE, n = 100)) > 0)
  expect_message(compute_crosstalk(c(1,3), g = g, use_ppi=FALSE))
})


