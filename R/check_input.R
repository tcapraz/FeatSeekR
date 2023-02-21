#' @title check_input
#'
#' @description Checks input data. Input data should be a 3 dimensional array
#' with samples x features x replicates or
#' list of replicate matrices or dataframes. Each matrix or dataframe
#' corresponds to a replicate and
#' should have matching dimensions (samples x features).
#' Column names should correspond to feature names
#' and must match between replicates.
#'
#' @param data input data provided to FeatSeek
#' @param replicates if data is a 2 dimensional array with samples x features
#' a dataframe with 2 columns named sample and replicates indicating which
#' sample corresponds to which replicate must be provided
#'
#'
#' @return input data reshaped to 3 dimensional array with
#' samples x features x replicates
#'
#' @keywords internal
check_input <- function(data, replicates, max_features){

    cnames <- colnames(data)
    if (all(vapply(cnames, function(x) is.null(x), logical(1)))) stop(
        "No feature names given or features not in correct dimension of data array!")

    if (!is.null(max_features) & max_features > dim(data)[2]) stop("Max features higher than features in data!")

    if( length(unique(replicates)) < 2) stop("At least 2 replicates required!")
    if (!length(replicates) == dim(data)[1]) stop("Replicate indicator vector not same length as samples in data!")
    if (any(table(replicates) < 2)) stop("Not every sample has at least 2 replicates!")
}

#' @title init_selected
#'
#' @description Checks if preselected init features are in input data.
#' If init is NULL, it is set to feature with highest replicate correlation.
#'
#' @param init preselected starting set of features
#' @param data input data
#' (3 dimensional array with samples x features x replicates)
#'
#'
#' @return names of initial set of feature
#'
#' @keywords internal
init_selected <- function(init, data, replicates){

    features <- dimnames(data)[[2]]
    if (is.null(features)) stop("No feature names given or features not in correct dimension of data array!")
    # check if init features are in data
    if(!is.null(init) & !all(init %in% features)){
        stop("Could not find init features in data!")
    }

    # start with feature with highest replicate correlation if no init
    # features were given
    if (is.null(init)){
        data <- data.frame(data)
        f <- vapply(seq_len(dim(data)[[2]]), function(x){
            m <- stats::lm(replicates~data[,x])
            s <- summary(m)
            s$fstatistic[1]
        }, numeric(1))
        names(f) <-  features
        init <- names(which.min(f))
    }
    init
}


#' @title filter
#'
#' @description This function removes features with mean pairwise replicate
#' pearson correlation < filter_thr
#'
#' @param data 3 dimensional array with samples x features x replicates.
#' @param filter_thr Mean pairwise replicate pearson correlation threshold
#'
#'
#' @return the filtered data array
#'
#' @keywords internal
filter <- function(data, init, filter_thr){
    message("Filtering out features with correlation < filter_thr
            across replicates!")
    # get mean of pairwise correlation between replicates,
    # i.e. the off-diagonal values of the correlation matrix
    cor <- apply(data,2, function(x)
        mean(stats::cor(x, use="pairwise.complete.obs")
             [upper.tri(!diag(nrow=dim(data)[3]))]))
    keep <-  cor > filter_thr
    features <- dimnames(data)[[2]]
    init_keep <- features %in% init
    features <- features[init_keep | keep]
    data <- data[, features,]
    p <- dim(data)[2] # get number of remaining features
    message("Remaining features: ", p)
    data
}

