test_that("fit lm returns data correctly",{
    data <-  array(rnorm(3*3*2), dim=c(3,3,2), dimnames=list(NULL, c("first", "second", "third"), NULL))
    v <- rnorm(3)
    m <- cbind(v,v*2, v*0.5)
    data_redundant <- abind::abind(m,m*0.1, along=3)

    S1 <- data[, 1, , drop = FALSE]
    # S1[, 1] = apply(data[, 1, , drop = FALSE], 1, mean, na.rm = TRUE)

    S2 <- data[, 1:2, , drop = FALSE]
    # S2[, 1] = apply(data[, 1, , drop = FALSE], 1, mean, na.rm = TRUE)
    # S2[, 2] = apply(data[, 2, , drop = FALSE], 1, mean, na.rm = TRUE)

    S3 <- data_redundant[, 1:2, , drop = FALSE]
    # S3[, 1] = apply(data_redundant[, 1, , drop = FALSE], 1, mean, na.rm = TRUE)
    # S3[, 2] = apply(data_redundant[, 2, , drop = FALSE], 1, mean, na.rm = TRUE)

    allzero <- array(0, dim=c(3,3,2), dimnames=list(NULL, c("v", "",  ""), NULL))

    # overwrite with zeros if features are collinear
    expect_identical(fit_lm(data_redundant, S3, 2), allzero)
    # return correct shape
    expect_identical(dim(fit_lm(data, S1, 2)), dim(data))
    expect_identical(dim(fit_lm(data, S2, 2)), dim(data))
    expect_identical(dim(fit_lm(data, S2, 3)), dim(data))
})
