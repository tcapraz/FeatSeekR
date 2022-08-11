#' @title variance_explained
#'
#' @param data 3 dimensional array samples x features x replicates
#' @param selected character vector of selected features
#'
#' @return average variance explained by selected features
#'
#' @keywords internal
variance_explained <- function(data,selected){
    r2 <- vapply(seq_len(dim(data)[[3]]), function(i){
        # calculate the fraction of variance explained by the selected feature set for each remaining feature.
        d0 <- data[,,i]
        rest <- d0
        data1 <- d0[,selected]
        # model remaining features by selected features
        model <- stats::lm(rest ~  data1 + 0)
        # model holds a fit for each remaining feature
        s <- suppressWarnings(summary(model))
        r <- mean(vapply(s, function(x){
          x$r.squared
        }, numeric(1)))
    }, numeric(1))
    # average across replicates
    mean(r2)
}

