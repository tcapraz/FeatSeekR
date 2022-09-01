test_that("Test calc_metric returns pairwise correlations", {
  data <-  array(c(1,2,3), dim=c(3,3,2), dimnames=list(NULL, c("first", "second", "third"), NULL))
  data[,1,1] <- c(3,2,1)
  data[,2,1] <- c(1,0.5,1)
  result <- c(-1,0,1)
  names(result) <- c("first", "second", "third")
  expect_equal(apply(data, 2, calcCor,  r = dim(data)[3]), result)

})

test_that("calcFstat returns correct output",{
  data <-  array(rnorm(100*3*2), dim=c(100,3,2), dimnames=list(NULL, c("first", "second", "third"), NULL))
  F_complete <- calcFstat(data, scale=TRUE, complete=TRUE)
  F_incomplete <- calcFstat(data, scale=TRUE, complete=FALSE)
  expect_identical(round(F_complete,5), round(F_incomplete,5))

})
