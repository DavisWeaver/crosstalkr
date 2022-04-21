
test_that("prep_stringdb returns an igraph", {
  skip("skipping due to long download times")
  #skip_if_offline()
  #skip_on_ci()
  expect_true(igraph::is.igraph(prep_stringdb()))
})
test_that("prep_biogrid returns an igraph", {
  skip("skipping due to long download times")
  #skip_if_offline()
  #skip_on_ci()
  expect_true(igraph::is.igraph(prep_biogrid()))
})


#have to skip on CI because CI can't find installed files
test_that("prep_stringdb works for multiple species", {
  expect_true(igraph::is.igraph(prep_stringdb(species = "pseudomonas aeruginosa")))
  expect_true(igraph::is.igraph(prep_stringdb(species = "burkholderia cepacia")))
  expect_true(igraph::is.igraph(prep_stringdb(species = "xanthomonas campestris")))
  expect_true(igraph::is.igraph(prep_stringdb(species = "pantoea agglomerans")))
})

test_that("prep_stringdb works when providing taxon id directly", {
  expect_true(igraph::is.igraph(prep_stringdb(species = "511145")))
  expect_true(igraph::is.igraph(prep_stringdb(species = 511145)))
})

test_that("prep_stringdb breaks when provided an unsupported species", {
  expect_error(prep_stringdb(species = "fluffy bunnies"))
  expect_error(prep_stringdb(species = "dinodinodino"))
  expect_error(prep_stringdb(species = "unicorn dragons"))

})

#test "to_taxon_id"
test_that("to_taxon_id works for a variety of species", {
  expect_equal(to_taxon_id("pseudomonas aeruginosa"), 287)
  expect_equal(to_taxon_id("PseudoMonas AeruGINOSA"), 287)
  expect_equal(to_taxon_id("homo sapiens"), 9606)
})

test_that("to_taxon_id gives error message for improperly names species", {
  expect_error(to_taxon_id("fluffy bunnies"))
  expect_error(to_taxon_id("dragons"))
  expect_error(to_taxon_id("unicorns"))
})

test_that("to_taxon_id gives error message for real species that aren't available", {
  expect_error(to_taxon_id("escherichia coli"))
})
