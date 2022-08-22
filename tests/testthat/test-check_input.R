
test_that("check_input returns correct array or appropriate error messages", {
    expected1 <- array(c(1,1,1,1,2,2,2,2), dim=c(2,2,2), dimnames=list(NULL, c("first", "second"), NULL))
    expected2 <- array(c(1,1,1,1,2,2,2,2,3,3,3,3), dim=c(2,2,3), dimnames=list(NULL, c("first", "second"), NULL))

    empty_list <- list()
    empty_mat <- matrix()
    mat1 <- matrix(1,2,2,dimnames=list(NULL, c("first", "second")))
    mat2 <- matrix(2,2,2,dimnames=list(NULL, c("first", "second")))
    mat3 <- matrix(3,2,2,dimnames=list(NULL, c("first", "second")))
    mat3_wrong <- matrix(3,2,2,dimnames=list(NULL, c("fir", "sec")))
    mat_list1 <- list(mat1, mat2)
    mat_list2 <- list(mat1, mat2, mat3)
    mat_list_wrong_names <- list(mat1,mat2,mat3_wrong)

    empty_df <- data.frame()
    df1 <- data.frame(first=c(1,1), second=c(1,1))
    df2 <- data.frame(first=c(2,2), second=c(2,2))
    df3 <- data.frame(first=c(3,3), second=c(3,3))
    df3_wrong <- data.frame(x=c(3,3), z=c(3,3))

    df_list1 <- list(df1, df2)
    df_list2 <-  list(df1, df2, df3)

    df_list_wrong_names <-  list(df1, df2, df3_wrong)
    mixed_list <- list(mat1, df1)

    empty_arr <- array(dim=c(2,2,2))
    arr <- array(1, dim=c(2,2,2), dimnames=list(NULL, c("first", "second"), NULL))
    arr_wrong_names <- array(1, dim=c(2,2,2), dimnames=list(NULL, NULL, c("first", "second")))
    arr_wrong_dims <- array(1, dim=c(2,2))


    data <- array(rnorm(2000), dim= c(100,20))
    colnames(data) <- paste("feature", seq_len(ncol(data)))
    replicates <- data.frame(sample=rep(seq_len(50),each=2), replicate=rep(c(1,2),50))
    data_expected1 <- array(data[seq(1,100,by=2),], dim = c(50,20,1))
    data_expected2 <- array(data[seq(2,100,by=2),], dim = c(50,20,1))
    data_expected <- abind::abind(data_expected1 ,data_expected2)
    expect_identical(check_input(data, replicates), data_expected)

    replicates <- data.frame(sample=rep(seq_len(50),2), replicate=rep(c(1,2),each=50))
    data_expected1 <- array(data[seq_len(50),], dim = c(50,20,1))
    data_expected2 <- array(data[seq(51,100),], dim = c(50,20,1))
    data_expected <- abind::abind(data_expected1 ,data_expected2)
    expect_identical(check_input(data, replicates), data_expected)


    data <- array(rnorm(150*20), dim= c(150,20))
    replicates <- data.frame(sample=rep(seq_len(50),3), replicate=rep(c(1,2,3),each=50))
    data_expected1 <- array(data[seq_len(50),], dim = c(50,20,1))
    data_expected2 <- array(data[seq(51,100),], dim = c(50,20,1))
    data_expected3 <- array(data[seq(101,150),], dim = c(50,20,1))
    data_expected <- abind::abind(data_expected1 ,data_expected2, data_expected3 )
    expect_identical(check_input(data, replicates), data_expected)

    replicates <- replicates[seq_len(148),]
    data <- data[seq_len(148),]
    data_expected1 <- array(data[seq_len(50),], dim = c(50,20,1))
    data_expected2 <- array(data[seq(51,100),], dim = c(50,20,1))
    data_expected3 <- array(dim = c(50,20,1))
    data_expected3[seq(1,48),,]  <- data[seq(101,148),]
    data_expected <- abind::abind(data_expected1 ,data_expected2, data_expected3 )
    same <- (check_input(data, replicates)  == data_expected) | (is.na(check_input(data, replicates)) & is.na(data_expected))
    expect_true(all(same))

    expect_identical(check_input(arr), arr)
    expect_identical(check_input(mat_list1), expected1)
    expect_identical(check_input(mat_list2), expected2)

    expect_identical(check_input(df_list1), expected1)
    expect_identical(check_input(df_list2), expected2)

    expect_error(check_input(mixed_list), "Not all elements in input list are either dataframes or matrices!")

    expect_error(check_input(mat_list_wrong_names), "Feature names are not the same for all replicates!")
    expect_error(check_input(df_list_wrong_names), "Feature names are not the same for all replicates!")

    expect_error(check_input(arr_wrong_dims), "Replicates not found. Please provide a dataframe indicating which sample corresponds to which replicate!")
    expect_error(check_input(arr_wrong_names), "No feature names given or features not in correct dimension of data array!")

    expect_error(check_input(empty_list), "Please provide a list of matrices/dataframes or a 3 dimensional array as input!")
    expect_error(check_input(empty_df), "Please provide a list of matrices/dataframes or a 3 dimensional array as input!")
    expect_error(check_input(empty_mat), "Replicates not found. Please provide a dataframe indicating which sample corresponds to which replicate!")


    expect_error(check_input("test"), "Please provide a list of matrices/dataframes or a 3 dimensional array as input!")
})

test_that("init returns feature with highest replicate correlation or issues warning if init features not found", {
    data <-  array(c(1,2,3), dim=c(3,3,2), dimnames=list(NULL, c("first", "second", "third"), NULL))
    data_no_names <- array(c(1,2,3), dim=c(3,3,2))
    data_wrong_names <-  array(c(1,2,3), dim=c(3,3,2), dimnames=list(NULL, c("fir", "sec", "th"), NULL))

    data[,1,1] <- c(3,2,1)
    data[,2,1] <- c(1,0.5,1)
    init <- "first"
    highest_cor <- "third"

    expect_error(init_selected(init, data_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(NULL, data_no_names), "No feature names given or features not in correct dimension of data array!")
    expect_error(init_selected(init, data_wrong_names), "Could not find init features in data!")

    expect_identical(init_selected(init, data), init)

    expect_identical(init_selected(init=NULL, data), highest_cor)
})

test_that("filter removes features with replicate correlation < filter_thr and keeps user defined init features" , {
    data <-  array(c(1,2,3), dim=c(3,3,2), dimnames=list(NULL, c("first", "second", "third"), NULL))
    data[,1,1] <- c(3,2,1)
    filter_thr <- 0.5
    init1 <- "first"
    init2 <-  "second"

    filtered_data <- data[,2:3,]

    expect_identical(filter(data, init1,filter_thr), data)
    expect_identical(filter(data, init2,filter_thr), filtered_data)
})


