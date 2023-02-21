
test_that("check_input returns appropriate error messages", {

    data <-  array(c(1,2,3), dim=c(6,3), dimnames=list(NULL, c("first", "second", "third")))
    data_no_names <- array(c(1,2,3), dim=c(6,3))
    data_wrong_names <-  array(c(1,2,3), dim=c(6,3), dimnames=list(NULL, c("fir", "sec", "th")))
    reps <- rep(c(1,2), each=3)
    wrong_reps <- rep(c(1,2), each=6)
    rep1 <- rep(1, 6)
    rep12 <- rep1
    rep12[1] <- 2
    expect_error(check_input(data_no_names, reps, max_features = 2), "No feature names given or features not in correct dimension of data array!")
    expect_error(check_input(data, wrong_reps, max_features = 2), "Replicate indicator vector not same length as samples in data!")
    expect_error(check_input(data, reps, max_features = 20), "Max features higher than features in data!")
    expect_error(check_input(data, rep1, max_features = 2), "At least 2 replicates required!")
    expect_error(check_input(data, rep12, max_features = 2), "Not every sample has at least 2 replicates!")


})

test_that("init returns feature with highest F-statistic or issues warning if init features not found", {
    data <-  array(c(1,2,3), dim=c(6,3), dimnames=list(NULL, c("first", "second", "third")))
    data_no_names <- array(c(1,2,3), dim=c(6,3))
    data_wrong_names <-  array(c(1,2,3), dim=c(6,3), dimnames=list(NULL, c("fir", "sec", "th")))
    reps <- rep(c(1,2), each=3)


    init <- "first"

    expect_error(init_selected(init, data_no_names,  reps), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(NULL, data_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(init, data_wrong_names), "Could not find init features in data!")

    expect_identical(init_selected(init, data), init)
    expect_identical(init_selected(init=NULL, data, reps), init)

})



