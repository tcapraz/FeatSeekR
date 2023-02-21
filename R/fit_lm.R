#' @title fit_lm
#'
#' @description Fit linear model for each feature as a function of
#' the previously selected features S.
#' The dimensions spanned by the selected features are
#' projected out of the data by setting each feature
#' to its residuals from the linear model fit.
#'
#' @param data data (3 dimensional array samples x features x replicates)
#' @param S set of selected features
#' @param k current iteration
#'
#' @return data with previously selected features projected out
#'
#' @keywords internal
fit_lm <- function(data, S, k){
    # get current features in data
    features <- dimnames(data)[[2]]
    d <- list()
    # fit linear model and overwrite data with residuals
    sel <- S[,seq_len(k-1)]
    model <- stats::lm(data ~ sel + 0, na.action = stats::na.exclude)
    eps <- array(NA, dim(data))
    eps[!is.na(data)] <- model$residuals
    data <- ifelse(abs(eps)<10^(-20), 0, eps)
    dimnames(data)[[2]] <- features
    data
}
