test_that("FeatSeek returns correct output", {
    data <-  array(rnorm(3*3*2), dim=c(3,3,2), dimnames=list(NULL, c("first", "second", "third"), NULL))
    v <- rnorm(3)
    m <- cbind(v,v*2, v*0.5)
    data_redundant <- abind::abind(m,m*0.1, along=3)
    dimnames(data_redundant) <- list(NULL, c("first", "second", "third"), NULL)
    init="first"
    max_features=2
    res1 <- FeatSeek(data, max_features)
    res2 <- FeatSeek(data, max_features, init)
    res3 <- FeatSeek(data_redundant, max_features, init)

    expect_identical(dim(res1)[1], as.integer(max_features))
    expect_true(inherits(res1,"data.frame"), "data.frame")
    expect_identical(res2$selected[[1]], init)

    # only one feature is selected as all the others are linear combination of init feature
    expect_identical(dim(res3)[1], as.integer(1))
})

