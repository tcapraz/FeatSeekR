
test_that("check_input returns appropriate error messages", {

    data <-  array(c(1,2,3), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))
    data_no_names <- array(c(1,2,3), dim=c(3,6))
    data_wrong_names <-  array(c(1,2,3), dim=c(3,6),
                dimnames=list(c("fir", "sec", "th"), NULL))
    reps <- rep(c(1,2), each=3)
    wrong_reps <- rep(c(1,2), each=6)
    rep1 <- rep(1, 6)
    rep12 <- rep1
    rep12[1] <- 2
    expect_error(check_input(data_no_names, reps, max_features = 2), "No feature names given or features not in correct dimension of data array!")
    expect_error(check_input(data, reps, max_features = 20), "Max features higher than features in data!")
    expect_error(check_input(data, rep1, max_features = 2), "At least 2 replicates required!")
    expect_error(check_input(data, rep12, max_features = 2), "Not every sample has at least 2 replicates!")


})

test_that("init returns feature with highest F-statistic or issues warning if init features not found", {
    data <-  array(rnorm(3*6), dim=c(3,6),
                dimnames=list(c("first", "second", "third"), NULL))
    data[1:2,1:3] <- 1
    data[1:2,4:6] <- 2
    data_no_names <- array(c(1,2,3), dim=c(3,6))
    data_wrong_names <-  array(c(1,2,3), dim=c(3,6),
                dimnames=list(c("fir", "sec", "th"), NULL))
    reps <- data.frame(replicates=rep(c(1,2), each=3))
    se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(data=data),
      colData=reps)
    se_no_names <- SummarizedExperiment::SummarizedExperiment(
      assays=list(data=data_no_names),
      colData=reps)
    se_wrong_names <- SummarizedExperiment::SummarizedExperiment(
      assays=list(data=data_wrong_names),
      colData=reps)

    init <- "first"
    expected_init <- "third"

    expect_error(init_selected(init, se_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(NULL, se_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(init, se_wrong_names), "Could not find init features in data!")

    expect_identical(init_selected(init, se), init)
    expect_identical(init_selected(init=NULL, se), expected_init)

})



