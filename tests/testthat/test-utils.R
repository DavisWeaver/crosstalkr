library(crosstalkr)
# Start by setting up test vectors of ensembl id
ensembl_test <- c("9606.ENSP00000000233", "9606.ENSP00000272298", "9606.ENSP00000371253",
                  "9606.ENSP00000253401")

#and entrez gene ids.
entrez_test <- c(79854, 148398, 26155, 84069, 57801, 9636,
                 54991, 254173, 8784)

test_that("as_gene_symbol converts ensemble ids to gene_symbol", {
  expect_equal(length(as_gene_symbol(ensembl_test)), length(ensembl_test))
  expect_equal(as_gene_symbol(ensembl_test), c("ARF5", "CALM2", "GART", "ARHGEF9"))
})

test_that("as_gene_symbol converts entrez ids to gene_symbol", {
  expect_equal(length(as_gene_symbol(entrez_test)), length(entrez_test))
  expect_equal(as_gene_symbol(entrez_test), c("LINC00115",
                                              "SAMD11",
                                              "NOC2L", "PLEKHN1", "HES4",
                                              "ISG15", "C1orf159", "TTLL10",
                                              "TNFRSF18" ))
})

test_that("is_ensembl returns true for ensembl ids and false otherwise", {
  expect_true(is_ensembl(ensembl_test))
  expect_true(is_ensembl(ensembl_test[1]))
  expect_true(is_ensembl("ENSG00000168374"))
  expect_true(is_ensembl("ENST00000168374"))
  expect_false(is_ensembl("ARF5"))
  expect_false(is_ensembl("entrez_test"))
})

test_that("is_entrez returns true for entrez_ids and false otherwise", {
  expect_true(is_entrez(entrez_test))
  expect_true(is_entrez(entrez_test[3]))
  expect_false(is_entrez("ARF5"))
  expect_false(is_entrez(ensembl_test))
})


test_that("entrez_type correctly identifies the type of entrez_id provided", {
  expect_equal(ensembl_type(ensembl_test), "PROTEINID")
  expect_equal(ensembl_type("9606.ENST00000000233"), "TRANSCRIPTID")
  expect_equal(ensembl_type("9606.ENSG00000000233"), "GENEID")
  expect_error(ensembl_type("ARF5"))
  expect_error(ensembl_type(1234))
})
