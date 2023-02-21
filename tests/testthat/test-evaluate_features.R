test_that("variance_explained returns correct output",{
    data <-  array(rnorm(50*30), dim=c(50,30), dimnames=list(NULL, paste("feature", seq_len(30))))
    selected <- paste("feature", seq_len(5))
    var_e <- variance_explained(data, selected)

    expect_true(inherits(var_e,"numeric"))
    expect_false(is.null(var_e))
})
