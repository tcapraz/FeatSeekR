
#
# test_that("reconstruction_error returns list of reconstruction errors", {
#     data <-  array(rnorm(50*30*2), dim=c(50,30,2), dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#
#     init="feature 1"
#
#     res1 <- FeatSeek(data, max_features=25)
#     res2 <- FeatSeek(data, max_features=1)
#
#     expect_identical(length(reconstruction_error(data, res1)), length(res1$selected))
#     expect_identical(length(reconstruction_error(data, res2)), length(res2$selected))
#     expect_true(inherits(reconstruction_error(data, res1),"list"))
#     expect_false(any(is.null(reconstruction_error(data, res1))))
# })


test_that("variance_explained returns correct output",{
    data <-  array(rnorm(50*30*2), dim=c(50,30,2), dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
    selected <- paste("feature", seq_len(5))
    var_e <- variance_explained(data, selected)

    expect_true(inherits(var_e,"numeric"))
    expect_false(is.null(var_e))
})
