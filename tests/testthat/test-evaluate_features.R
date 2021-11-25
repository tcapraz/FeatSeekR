test_that("get_opt returns optimal number of features",{
    data1 <-  array(rnorm(20*10*2), dim=c(20,10,2), dimnames=list(NULL, paste("feature", seq_len(10)), NULL))
    data2 <-  array(rnorm(50*30*2), dim=c(50,30,2), dimnames=list(NULL, paste("feature", seq_len(30)), NULL))

    init="feature 1"

    res1 <- FeatSeek(data1, max_features=2)
    res2 <- FeatSeek(data2, max_features=25)

    expect_message(get_opt(res1), "Less than 4 features selected! Cannot determine optimal set, returning number of selected features!")
    expect_identical(get_opt(res1), as.integer(2))
    expect_identical(typeof(get_opt(res2)), "integer")
})


test_that("svd_entropy returns list of svd entropies", {
    data <-  array(rnorm(50*30*2), dim=c(50,30,2), dimnames=list(NULL, paste("feature", seq_len(30)), NULL))

    init="feature 1"

    res1 <- FeatSeek(data, max_features=25)
    res2 <- FeatSeek(data, max_features=1)

    expect_identical(length(svd_entropy(data, res1)), length(res1$selected))
    expect_identical(length(svd_entropy(data, res2)), length(res2$selected))
    expect_true(inherits(svd_entropy(data, res1),"list"))
    expect_false(any(is.null(svd_entropy(data, res1))))

})

test_that("reconstruction_error returns list of reconstruction errors", {
    data <-  array(rnorm(50*30*2), dim=c(50,30,2), dimnames=list(NULL, paste("feature", seq_len(30)), NULL))

    init="feature 1"

    res1 <- FeatSeek(data, max_features=25)
    res2 <- FeatSeek(data, max_features=1)

    expect_identical(length(reconstruction_error(data, res1)), length(res1$selected))
    expect_identical(length(reconstruction_error(data, res2)), length(res2$selected))
    expect_true(inherits(reconstruction_error(data, res1),"list"))
    expect_false(any(is.null(reconstruction_error(data, res1))))
})
