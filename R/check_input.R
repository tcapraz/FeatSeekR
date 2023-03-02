#' @title check_input
#'
#' @description Checks input data. Input data should be a 2 dimensional array
#' with samples x features or SummarizedExperiment carrying one assay named data
#' and colData indicating which sample belongs to which replicate
#'
#'
#' @param data input data provided to FeatSeek
#' @param replicates if data is a 2 dimensional array with samples x features
#' a vector indicating which sample corresponds to which replicate
#' must be provided
#'
#'
#' @return input data reshaped to 3 dimensional array with
#' samples x features x replicates
#'
#' @keywords internal
check_input <- function(data, max_features, replicates=NULL){
    if (!inherits(data,"SummarizedExperiment")){
        reps <- data.frame(replicates =replicates)
        se <- SummarizedExperiment::SummarizedExperiment(assays=list(data=data), colData=reps)
    } else {
        se <- data
    }
    fnames <- rownames(assays(se)$data)
    if (all(vapply(fnames, function(x) is.null(x), logical(1)))) stop(
        "No feature names given or features not in correct dimension of data array!")

    if (!is.null(max_features) & max_features > dim(assays(se)$data)[1]) stop("Max features higher than features in data!")

    if( length(unique(colData(se)$replicates)) < 2) stop("At least 2 replicates required!")
    if (!length(colData(se)$replicates) == dim(assays(se)$data)[2]) stop("Replicate indicator vector not same length as samples in data!")
    if (any(table(colData(se)$replicates) < 2)) stop("Not every sample has at least 2 replicates!")
    se
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
init_selected <- function(init, se){

    features <- dimnames(assays(se)$data)[[1]]
    replicates <- colData(se)$replicates
    if (is.null(features)) stop("No feature names given or features not in correct dimension of data array!")
    # check if init features are in data
    if(!is.null(init) & !all(init %in% features)){
        stop("Could not find init features in data!")
    }

    # start with feature with highest replicate correlation if no init
    # features were given
    if (is.null(init)){
        data <- t(data.frame(assays(se)$data))
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

