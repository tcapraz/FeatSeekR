

test_that("calcFstat returns correct output",{
    data <-  array(rnorm(100*3), dim=c(100,3), dimnames=list(NULL, c("first", "second", "third")))
    conds <- rep(seq_len(50), 2)
    
    expect_error(calcFstat(data, conds))
    
    conds <- as.factor(conds)
    Fstat <- calcFstat(data, conds)
    
    expect_error(calcFstat(data, c(1,2)))

    expect_true(is(Fstat, "numeric"))
    expect_identical(length(Fstat), dim(data)[[2]])
})
