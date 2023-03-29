test_that("fit lm returns data correctly",{
    data <-  array(rnorm(6*3), dim=c(6,3), dimnames=list(NULL, c("first", "second", "third")))
    v <- rnorm(6)
    data_redundant <- cbind(v,v*2, v*0.5)

    S1 <- data[, 1,  drop = FALSE]
    # S1[, 1] = apply(data[, 1, , drop = FALSE], 1, mean, na.rm = TRUE)

    S2 <- data[, 1:2,  drop = FALSE]
    # S2[, 1] = apply(data[, 1, , drop = FALSE], 1, mean, na.rm = TRUE)
    # S2[, 2] = apply(data[, 2, , drop = FALSE], 1, mean, na.rm = TRUE)

    S3 <- data_redundant[, 1:2,  drop = FALSE]
    # S3[, 1] = apply(data_redundant[, 1, , drop = FALSE], 1, mean, na.rm = TRUE)
    # S3[, 2] = apply(data_redundant[, 2, , drop = FALSE], 1, mean, na.rm = TRUE)
    out <- fit_lm(data, S1, 2)
    
    
    expect_error(fit_lm(data, 1, 2))
    expect_error(fit_lm(data, S1, c(1,2)))
    expect_error(fit_lm(data, S1, 10))
    
    expect_true(is(out, "array"))
    expect_true(all(colnames(out) == c("first", "second", "third")))
    
    # return correct shape
    expect_identical(dim(fit_lm(data, S1, 2)), dim(data))
    expect_identical(dim(fit_lm(data, S2, 2)), dim(data))
    expect_identical(dim(fit_lm(data, S2, 3)), dim(data))
})
