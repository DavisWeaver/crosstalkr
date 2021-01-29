
library(crosstalkr)
# Start by testing if ensembl_convert works
ensembl_test <- c("9606.ENSP00000000233", "9606.ENSP00000272298", "9606.ENSP00000371253",
                  "9606.ENSP00000253401")

test_that("ensembl_convert works", {
  expect_equal(length(ensembl_convert(ensembl_test)), length(ensembl_test))
  expect_equal(ensembl_convert(ensembl_test), c("ARF5", "CALM2", "GART", "ARHGEF9"))
})
ensembl_convert(ensembl_test)

test_that("prep_stringdb returns an igraph", {
  expect_true(igraph::is.igraph(prep_stringdb()))
})


test_that("prep_biogrid returns an igraph", {
  expect_true(igraph::is.igraph(prep_biogrid()))
})
