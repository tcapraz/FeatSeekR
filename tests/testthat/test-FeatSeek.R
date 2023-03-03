test_that("FeatSeek returns correct output", {
    data <-  array(rnorm(6*3), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))

    init="first"
    reps <- rep(c(1,2), each=3)
    max_features=2
    res1 <- FeatSeek(data, reps, max_features=max_features)
    res2 <- FeatSeek(data, reps, max_features=max_features, init=init)

    expect_identical(dim(res1)[1], as.integer(max_features))
    expect_true(inherits(res1,"SummarizedExperiment"), "SummarizedExperiment")
    expect_identical(SummarizedExperiment::rowData(res2)$selected[[1]], init)

})

