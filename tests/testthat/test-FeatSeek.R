test_that("FeatSeek returns correct output", {
    data <-  array(rnorm(6*3), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))
    
    data_wrong <- array(rnorm(6*3), dim=c(3,6))
    
    init="first"
    conds <- rep(seq_len(3), 2)
    max_features=2
    res1 <- FeatSeek(data, conds, max_features=max_features)
    res2 <- FeatSeek(data, conds, max_features=max_features, init=init)
    
    conds_wrong <- seq_len(10)
    
    char_conds <- c("first", "first", "second", "second")
    num_init <- 1
    max_features_wrong <- "ten" 
      
    expect_error(FeatSeek(data, char_conds, max_features=2))
    expect_error(FeatSeek(data, conds, init=num_init, max_features=2))
    expect_error(FeatSeek(data, conds, max_features=max_features_wrong))
    
    expect_error(FeatSeek(data, conds, max_features=10), "Max features higher than features in data!")
    expect_error(FeatSeek(data, conds_wrong, max_features=2), "Condition factor not same length as samples in data!")
    expect_error(FeatSeek(data_wrong, conds, max_features=2), "No feature names given or features not in correct dimension of data array!")
    
    expect_identical(dim(res1)[1], as.integer(max_features))
    expect_true(is(res1,"SummarizedExperiment"), "SummarizedExperiment")
    expect_identical(SummarizedExperiment::rowData(res2)$selected[[1]], init)
    expect_true(length(SummarizedExperiment::rowData(res2)$selected) == max_features)
})

