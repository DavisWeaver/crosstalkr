
test_that("prep_stringdb returns an igraph", {
  skip_on_ci()
  expect_true(igraph::is.igraph(prep_stringdb()))
})
test_that("prep_biogrid returns an igraph", {
  skip_on_ci()
  expect_true(igraph::is.igraph(prep_biogrid()))
})

