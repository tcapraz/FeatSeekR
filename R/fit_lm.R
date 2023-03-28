#' @title fit_lm
#'
#' @description Fit linear model for each feature as a function of
#' the previously selected features S.
#' The dimensions spanned by the selected features are
#' projected out of the data by setting each feature
#' to its residuals from the linear model fit.
#'
#' @param data \code{data} (2 dimensional array samples x features)
#' @param S set of selected features
#' @param k current iteration
#'
#' @return \code{data} with previously selected features projected out
#' 
#' @importFrom stats lm na.exclude
#' @importFrom methods is
#' @keywords internal
fit_lm <- function(data, S, k){
    stopifnot(
        is.numeric(data),
        is.numeric(S),
        k-1 <= ncol(S),
        length(k) == 1
    )
    # get current features in data
    features <- dimnames(data)[[2]]
    d <- list()
    # fit linear model and overwrite data with residuals
    sel <- S[, seq_len(k-1)]
    model <- lm(data ~ sel + 0, na.action=na.exclude)
    eps <- array(NA, dim(data))
    eps[!is.na(data)] <- model$residuals
    data <- ifelse(abs(eps) < 10^(-20), 0, eps)
    dimnames(data)[[2]] <- features
    data
}
