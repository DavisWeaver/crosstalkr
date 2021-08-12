
if(curl::has_internet()) {
  test_that("prep_stringdb returns an igraph", {
    expect_true(igraph::is.igraph(prep_stringdb()))
  })


  test_that("prep_biogrid returns an igraph", {
    expect_true(igraph::is.igraph(prep_biogrid()))
  })

}
