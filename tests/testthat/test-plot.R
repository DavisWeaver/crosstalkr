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

df_true <- compute_crosstalk(c(1,3), g = g, use_ppi=FALSE)
g <- igraph::sample_gnp(n = 1000, p = 500/1000)

# test_that("check_crosstalk identifies crosstalk dfs as true", {
#   expect_true(check_crosstalk(df_true))
#   expect_true(check_crosstalk(compute_crosstalk(c(1), g = g, use_ppi=FALSE, n = 100)))
# })

test_that("check_crosstalk identifies returns false for incorrect input", {
  expect_false(check_crosstalk(1))
  expect_false(check_crosstalk("some string"))

})
