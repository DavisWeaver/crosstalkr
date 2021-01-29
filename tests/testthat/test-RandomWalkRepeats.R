library(crosstalkr)

v1 = (c(1,1,1,0))
v2 = c(0,0,0,1)
v3 = c(1,1,1,0)
v4 = c(0,0,0,1)

w = matrix(data = c(v1,v2,v3,v4), ncol = 4, nrow = 4)

seed = 1
p <- sparseRWR(w, seed)

test_that("sparseRWR returns a list of 3 ouputs of the correct types and lengths", {
  expect_true(is.double(p[[1]]))
  expect_equal(length(p[[1]]), 4)
  expect_true(is.double(p[[2]]))
  expect_equal(length(p[[2]]), 1)
  expect_equal(length(p[[3]]), 1)
  expect_true(is.integer(p[[3]]))
})

test_that("sparseRWR completes in < 1000 iterations", {
  expect_lt(p[[3]], 1000)
})

test_that("sparseRWR output equals expected", {
  expect_equal(sparseRWR(w, 1)[[1]], c(0.709, 0, 0.109, 0), tolerance = 0.001)
  expect_equal(sparseRWR(w, 2)[[1]], c(0.1090909, 0.6, 0.1090909, 0), tolerance = 0.001)
  expect_equal(sparseRWR(w, 3)[[1]], c(0.1090909, 0, 0.7090909, 0), tolerance = 0.001)
  expect_equal(sparseRWR(w, 4)[[1]], c(0.07272, 0.4, 0.07272, 1), tolerance = 0.001)
  expect_equal(sparseRWR(w, c(1,2))[[1]], c(0.40909, 0.3, 0.10909, 0), tolerance = 0.001)
  expect_lte(sum(sparseRWR(w, c(1,3))[[1]]), 1)
})

#Now need to make sure that it can be run in parallel
#Previously we had a bug where it would fly an error for the first set of workers.
library(foreach)

m <- replicate(1000, sample(x = c(0,1), size = 1000, replace = TRUE))
w <- Matrix::Matrix(m, sparse = TRUE)
w <- Matrix::t(Matrix::t(w)/Matrix::colSums(w)) #normalize based on the column sum.
seeds <- sample(1:nrow(w), size = 32)

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

n = 8
null_dist <-
  foreach::foreach(i = 1:n, .errorhandling = 'pass', .packages = "Matrix") %dopar% {
    crosstalkr::sparseRWR(w, seed_proteins = seeds, norm = FALSE)[[1]]
  }

test_that("parallel execution doesn't return an error", {
  expect_true(is.numeric(null_dist[[1]]))
})
