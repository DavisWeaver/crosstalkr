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
})

g <- igraph::sample_gnp(n = 1000, p = 300/1000)

test_that("compute_crosstalk identifies crosstalkers for larger graphs", {
  expect_true(is.data.frame(compute_crosstalk(c(1,3,5,8,10), g = g, use_ppi = FALSE, n = 100)))
  expect_true(nrow(compute_crosstalk(c(1,3,5,8,10), g = g, use_ppi = FALSE, n = 100)) > 0)
  expect_message(compute_crosstalk(c(1,3), g = g, use_ppi=FALSE))
})

test_that("compute_crosstalk doesn't break when you provide invalid vertez names", {
  expect_true(is.data.frame(compute_crosstalk(c(1,3,5,8,10, 1002), g = g, use_ppi = FALSE, n = 100)))
})


test_that("compute crosstalk works for non-human vertices", {
  skip_if_offline()
  expect_true(is.data.frame(compute_crosstalk("511145.b0003", species = 511145)))
})

test_that("gfilter.ct returns the expected subnetwork", {
  g <- gfilter.ct(seeds=c(1,3,5,8,10), g=g, use_ppi=FALSE, n =100)
  expect_true(length(igraph::V(g)) < 1000)
  expect_true(length(igraph::V(g)) >= 5) #at least the seeds should have come up as significant
})

g <- igraph::sample_gnp(n = 1000, p = 300/1000)
val <- sample.int(500, size = 1000, replace = TRUE)
names(val) <- 1:1000
test_that("gfilter.value returns the expected number of vertices", {
  obj = gfilter.value(g=g, val=val, n = 100, use_ppi = FALSE, desc = TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  obj = gfilter.value(g=g, val = 1:1000, n=100, use_ppi = FALSE, desc=TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "value")
  expect_true(all(val >= 900))
})

test_that("gfilter correctly calls gfilter.value", {
  obj = gfilter(g=g, method = "value", val=val, n = 100, use_ppi = FALSE, desc = TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  obj = gfilter(g=g, method = "value", val = 1:1000, n=100, use_ppi = FALSE, desc=TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "value")
  expect_true(all(val >= 900))
})

nps <- abs(calc_np_all(exp=val, g = g))
nps <- sort(nps, decreasing=TRUE)
min <- min(nps[1:100])
test_that("gfilter.np returns the expected number of vertices", {
  obj = gfilter.np(g=g,  use_ppi = FALSE, val=val, n = 100, desc = TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  obj = gfilter.np(g=g, val = 1:1000, n=100, use_ppi = FALSE, desc=TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "np")
  expect_true(all(val >= min))
})

test_that("gfilter correctly calls gfilter.np", {
  obj = gfilter(g=g, method = "np",  use_ppi = FALSE, val=val, n = 100, desc = TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  obj = gfilter(g=g, method = "np", val = 1:1000, n=100, use_ppi = FALSE, desc=TRUE)
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "np")
  expect_true(all(val >= min))
})

####Now we'll repeat the trick above for various node ranking igraph methods###

##Starting with igraph::degree
degree <- igraph::degree(g)
degree <- sort(degree, decreasing=TRUE)
min <- min(degree[1:100])
test_that("gfilter.igraph_method works for degree", {
  obj = gfilter.igraph_method(g=g,  use_ppi = FALSE, method=igraph::degree,
                              n = 100, desc = TRUE, val_name = "degree")
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "degree")
  expect_true(all(val >= min))
})

test_that("gfilter correctly calls gfilter.igraph_method for degree", {
  obj = gfilter(g=g, use_ppi = FALSE, igraph_method=igraph::degree,
                n = 100, desc = TRUE, val_name = "degree")
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "degree")
  expect_true(all(val >= min))
})

##now igraph::betweenness
between <- igraph::betweenness(g)
between <- rev(sort(between))
min <- min(between[1:100])
test_that("gfilter.igraph_method works for betweenness", {
  obj = gfilter.igraph_method(g=g,  use_ppi = FALSE, method=igraph::betweenness,
                              n = 100, desc = TRUE, val_name = "betweenness")
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "betweenness")
  expect_true(all(val >= min))
})

test_that("gfilter correctly calls gfilter.igraph_method for betweenness", {
  obj = gfilter(g=g, use_ppi = FALSE, igraph_method=igraph::betweenness,
                n = 100, desc = TRUE, val_name = "degree")
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "degree")
  expect_true(all(val >= min))
})

test_that("gfilter correctly calls gfilter.igraph_method for betweenness (specified as a character", {
  obj = gfilter(g=g, use_ppi = FALSE, igraph_method="betweenness",
                n = 100, desc = TRUE, val_name = "degree")
  expect_equal(length(igraph::V(obj)), 100)
  val = igraph::get.vertex.attribute(obj, name = "degree")
  expect_true(all(val >= min))
})

##now igraph::
