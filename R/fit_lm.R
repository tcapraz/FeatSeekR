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
    r <- dim(data)[[3]]
    d <- list()
    # fit linear model for each replicate and overwrite data with residuals
    sel <- S[,seq_len(k-1), drop = FALSE]
    for (j in seq_len(r)){
        model = stats::lm(data[, , j] ~ sel + 0, na.action = stats::na.exclude)
        eps <- array(NA, dim(data[, , j]))
        eps[!is.na(data[, , j])] <- model$residuals
        data[, , j] = ifelse(abs(eps)<10^(-10), 0, eps)
    }
    data
}
