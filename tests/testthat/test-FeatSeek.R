test_that("FeatSeek returns correct output", {
    data <-  array(rnorm(6*3), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))
    
    data_wrong <- array(rnorm(6*3), dim=c(3,6))
    
    init="first"
    reps <- rep(c(1,2), each=3)
    max_features=2
    res1 <- FeatSeek(data, reps, max_features=max_features)
    res2 <- FeatSeek(data, reps, max_features=max_features, init=init)
    
    reps_wrong <- seq_len(10)
    
    char_reps <- c("first", "first", "second", "second")
    num_init <- 1
    max_features_wrong <- "ten" 
      
    expect_error(FeatSeek(data, char_reps, max_features=2))
    expect_error(FeatSeek(data, reps, init=num_init, max_features=2))
    expect_error(FeatSeek(data, reps, max_features=max_features_wrong))
    
    expect_error(FeatSeek(data, reps, max_features=10), "Max features higher than features in data!")
    expect_error(FeatSeek(data, reps_wrong, max_features=2), "Replicate indicator vector not same length as samples in data!")
    expect_error(FeatSeek(data_wrong, reps, max_features=2), "No feature names given or features not in correct dimension of data array!")
    
    expect_identical(dim(res1)[1], as.integer(max_features))
    expect_true(is(res1,"SummarizedExperiment"), "SummarizedExperiment")
    expect_identical(SummarizedExperiment::rowData(res2)$selected[[1]], init)
    expect_true(length(SummarizedExperiment::rowData(res2)$selected) == max_features)
})

