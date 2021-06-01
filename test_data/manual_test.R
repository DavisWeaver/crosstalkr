
#lets benchmark cpp_RWR  compared to sparse_RWR
g <- prep_stringdb(cache = "test_data", min_score = 600)
load("./test_data/seed_proteins.Rda")
w <- igraph::as_adjacency_matrix(g)
rownames <- rownames(w)
w <- norm_colsum(w)

microbenchmark::microbenchmark(cpp_RWR(w =w, rownames = rownames, seed_proteins = seed_proteins),
                               sparseRWR(seed_proteins = seed_proteins, w = w))

g2 <- prep_stringdb(cache = NULL, min_score = 400)
#Test with multiple cores
time1 <- Sys.time()
test1 <- bootstrap_null(g = g, ncores = 4, seed_proteins = seed_proteins,
                       cache = "test_data", n = 10000)

time2 <- Sys.time()
diff1 <- time2-time1


#Test with one core.l
time1 <- Sys.time()
profvis::profvis(expr = bootstrap_null(g = g, ncores = 1, seed_proteins = seed_proteins,
                                       cache = "test_data", n = 10000))

time2 <- Sys.time()
diff2 <- time2-time1

#test compute_crosstalk (high-level function)

test3 <- compute_crosstalk(seed_proteins = seed_proteins, ppi = "stringdb",
                           cache = "test_data", seed_name = "angiogenesis",
                           ncores = 6)

#see if the graph function works

plot_ct(test3, g = g, prop_keep = 0.2, label_prop = 0.2)


#test whether we can actually run biogrid - we cannot

test4 <- prep_biogrid(cache = "test_data")
w <- igraph::as_adjacency_matrix(test4)
test5 <- sparseRWR(w, seed_proteins = seed_proteins, norm= FALSE)

#see if we can run sparseRWR on random seeds.
w <- igraph::as_adjacency_matrix(g)
test6 <- sparseRWR(w, seed_proteins = random_seeds)

#save outfile
my_log <- file("./test_data/my_log.txt") # File name of output log

sink(my_log, append = TRUE, type = "output") # Writing console output to log file
sink(my_log, append = TRUE, type = "message")
closeAllConnections() # Close connection to log file


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
