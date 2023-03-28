#' @title check_input
#'
#' @description Checks input \code{data}. Input \code{data} should be a 2 dimensional \code{array}
#' with samples x features or \code{SummarizedExperiment} carrying one assay named \code{data}
#' and \code{colData} indicating which sample belongs to which replicate
#'
#'
#' @param data input \code{data} provided to \code{FeatSeek} either \code{SummarizedExperiment} or
#' 2 dimensional \code{array} with features x samples
#' @param replicates if \code{data} is a 2 dimensional \code{array} with features x samples
#' a vector indicating which sample corresponds to which replicate
#' must be provided
#'
#'
#' @return \code{SummarizedExperiment} where replicate information is stored in colData
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @keywords internal
check_input <- function(data, max_features, replicates=NULL){
    if (!inherits(data,"SummarizedExperiment")){
        reps <- data.frame(replicates=replicates)
        se <- SummarizedExperiment(assays=list(data=data), colData=reps)
    } else {
        se <- data
    }
    fnames <- rownames(se)
    if (all(vapply(fnames, function(x) is.null(x), logical(1)))) stop(
        "No feature names given or features not in correct dimension of data array!")

    if (!is.null(max_features) & max_features > nrow(se)) stop("Max features higher than features in data!")

    if( length(unique(se$replicates)) < 2) stop("At least 2 replicates required!")
    if (!length(se$replicates) == ncol(se)) stop("Replicate indicator vector not same length as samples in data!")
    if (any(table(se$replicates) < 2)) stop("Not every sample has at least 2 replicates!")
    se
}

#' @title init_selected
#'
#' @description Checks if preselected init features are in input data.
#' If init is NULL, it is set to feature with highest replicate correlation.
#'
#' @param init preselected starting set of features
#' @param data input data as SummarizedExperiment
#'
#'
#'
#' @return names of initial set of feature
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom stats lm
#' @importFrom methods is
#' 
#' @keywords internal
init_selected <- function(init, se){
    stopifnot(
        is(se, "SummarizedExperiment"),
        # 'n_features' should be a scalar integer and 
        # mustn't exeed the number of available features
        is.character(init) | is.null(init)
    )
    features <- dimnames(assay(se, "data"))[[1]]
    replicates <- se$replicates
    if (is.null(features)) stop("No feature names given or features not in correct dimension of data array!")
    # check if init features are in data
    if(!is.null(init) & !all(init %in% features)){
        stop("Could not find init features in data!")
    }

    # start with feature with highest replicate correlation if no init
    # features were given
    if (is.null(init)){
        data <- t(data.frame(assay(se, "data")))
        f <- vapply(seq_len(dim(data)[[2]]), function(x){
            m <- lm(replicates ~ data[,x])
            s <- withCallingHandlers(summary(m),
                                warning = function(w){
                                    if(startsWith(conditionMessage(w),
                                        "essentially perfect fit")){
                                        invokeRestart("muffleWarning")
                                    } else {
                                        message(w$message)
                                    }
                                })
            s$fstatistic[1]
        }, numeric(1))
        names(f) <-  features
        init <- names(which.min(f))
    }
    init
}



