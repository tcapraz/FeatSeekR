#' @title check_input
#'
#' @description Checks input \code{data}. 
#' Input \code{data} should be a 2 dimensional \code{array}
#' with features x samples or \code{SummarizedExperiment} carrying one assay 
#' named \code{data} and \code{colData} indicating which sample belongs 
#' to which condition
#'
#'
#' @param data input \code{data} provided to \code{FeatSeek} either 
#' \code{SummarizedExperiment} or
#' 2 dimensional \code{array} with features x samples
#' @param conditions if \code{data} is a 2 dimensional \code{array} 
#' with features x samples
#' a factor indicating which sample corresponds to which condition
#' must be provided
#'
#'
#' @return \code{SummarizedExperiment} where condition information is stored in 
#' colData
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @keywords internal
check_input <- function(data, max_features, conditions=NULL){
    if (!is(data,"SummarizedExperiment")){
        cond<- data.frame(conditions=as.factor(conditions))
        if (!length(cond$conditions) == ncol(data)) 
            stop(strwrap(
                "Condition factor not same length as samples in 
                data!")
            )
        se <- SummarizedExperiment(assays=list(data=data), colData=cond)
    } else {
        se <- data
    }
    fnames <- rownames(se)
    if (all(vapply(fnames, function(x) is.null(x), logical(1)))) 
        stop(strwrap(prefix = " ", initial = "",
            "No feature names given or features not in correct dimension of 
            data array!")
        )
    if (!is.null(max_features) & max_features > nrow(se)) 
        stop("Max features higher than features in data!")
    if( length(unique(se$conditions)) < 2) 
        stop("At least 2 conditions required!")
    if (any(table(se$conditions) < 2)) 
        stop("Not every condition has at least 2 replicates!")
    se
}

#' @title init_selected
#'
#' @description Checks if preselected init features are in input data.
#' If init is NULL, it is set to feature with highest condition correlation.
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
    conditions <- se$conditions
    if (is.null(features)) 
        stop(strwrap(prefix = " ", initial = "",
            "No feature names given or features not in correct dimension of 
            data array!")
        )
    # check if init features are in data
    if(!is.null(init) & !all(init %in% features)){
        stop("Could not find init features in data!")
    }

    # start with feature with highest condition correlation if no init
    # features were given
    if (is.null(init)){
        data <- t(assay(se, "data"))
        f <- calcFstat(data, se$conditions)
        names(f) <-  features
        init <- names(which.max(f))
    }
    init
}



