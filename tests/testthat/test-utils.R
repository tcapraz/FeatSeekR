

test_that("calcFstat returns correct output",{
    data <-  array(rnorm(100*3), dim=c(100,3), dimnames=list(NULL, c("first", "second", "third")))
    reps <- rep(c(1,2), each=50)
    Fstat <- calcFstat(data, reps, scale=TRUE)
    expect_identical(length(Fstat), dim(data)[[2]])

})
