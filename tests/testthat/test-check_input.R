
test_that("check_input returns appropriate error messages", {

    data <-  array(c(1,2,3), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))
    data_no_names <- array(c(1,2,3), dim=c(3,6))
    data_wrong_names <-  array(c(1,2,3), dim=c(3,6),
                dimnames=list(c("fir", "sec", "th"), NULL))
    conds <- rep(seq_len(3), 2)
    wrong_conds <- rep(seq_len(6), 2)
    cond1 <- rep(1, 6)
    cond12 <- cond1
    cond12[1] <- 2
    long_conds <- rep(c(1,2), each=10)
    
    expect_error(check_input(data_no_names, conds, max_features = 2), "No feature names given or features not in correct dimension of data array!")
    expect_error(check_input(data, conds, max_features = 20), "Max features higher than features in data!")
    expect_error(check_input(data, cond1, max_features = 2), "At least 2 conditions required!")
    expect_error(check_input(data, cond12, max_features = 2), "Not every condition has at least 2 replicates!")
    expect_error(check_input(data, long_conds, max_features=2), "Condition factor not same length as samples in data!")
    expect_true(is(check_input(data, conds, max_features = 2), "SummarizedExperiment"))
    
    se <- check_input(data, conds, max_features = 2)
    expect_true(all(se$conditions == conds))
    expect_true(all(rownames(se) == c("first", "second", "third")))

})

test_that("init returns feature with highest F-statistic or issues warning if init features not found", {
    data <-  array(rnorm(3*6), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))
    data[1:2,1:3] <- 1
    data[1:2,4:6] <- 2
    data_no_names <- array(c(1,2,3), dim=c(3,6))
    data_wrong_names <-  array(c(1,2,3), dim=c(3,6),
                dimnames=list(c("fir", "sec", "th"), NULL))
    conds <- data.frame(conditions=as.factor(rep(seq_len(3), 2)))
    se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(data=data),
      colData=conds)
    se_no_names <- SummarizedExperiment::SummarizedExperiment(
      assays=list(data=data_no_names),
      colData=conds)
    se_wrong_names <- SummarizedExperiment::SummarizedExperiment(
      assays=list(data=data_wrong_names),
      colData=conds)

    init <- "first"
    expected_init <- "third"

    expect_error(init_selected(init, se_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(NULL, se_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(init, se_wrong_names), "Could not find init features in data!")
    
    expect_true(is(init_selected(init, se), "character"))
    expect_true(length(init_selected(init, se)) == 1)
    expect_identical(init_selected(init, se), init)
    expect_identical(init_selected(init=NULL, se), expected_init)

})



