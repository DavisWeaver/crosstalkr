##Test the combine_Null function
load(system.file("test_data/toy_combine.Rda", package = "disruptr"))

#first case when all the provided arguments are un-aggregated
arg1 <- df
arg2 <- df
arg3 <- df
arg0 <- combine_null(arg1,arg2,arg3)
test_that("combine_null aggregates multiple inner loop tasks"{
  expect_equal(ncol(arg0), 5)
  expect_true(all(arg0$n == 3))
})


#second case when there is a provided argument that is aggregated first
test2 <- combine_null(arg0,arg1,arg2,arg3)

test_that("combine_null aggregates multiple inner loop tasks and one pre-aggregated data.frame"{
  expect_equal(ncol(test2), 5)
  expect_true(all(test2$n == 6))
})
