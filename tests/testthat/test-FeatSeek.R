test_that("FeatSeek returns correct output", {
    data <-  array(rnorm(6*3), dim=c(6,3), dimnames=list(NULL, c("first", "second", "third")))

    init="first"
    reps <- rep(c(1,2), each=3)
    max_features=2
    res1 <- FeatSeek(data, reps, max_features=max_features)
    res2 <- FeatSeek(data, reps, max_features=max_features, init=init)

    expect_identical(dim(res1)[1], as.integer(max_features))
    expect_true(inherits(res1,"data.frame"), "data.frame")
    expect_identical(res2$selected[[1]], init)

})

