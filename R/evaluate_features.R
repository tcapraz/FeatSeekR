#' @title variance_explained
#'
#' @param data 3 dimensional array samples x features x replicates
#' @param selected character vector of selected features
#'
#' @return average variance explained by selected features
#'
#' @keywords internal
variance_explained <- function(data,selected){

    # calculate the fraction of variance explained by the selected feature set for each remaining feature.
    rest <- data
    data_sel <- data[,selected]
    # model remaining features by selected features
    model <- stats::lm(rest ~  data_sel)
    # model holds a fit for each remaining feature
    # catch warning that fit is perfect, as this is expected at certain number
    # of selected features
    s <- withCallingHandlers(summary(model),
                            warning = function(w){
                                if(startsWith(conditionMessage(w), "essentially perfect fit")){
                                    invokeRestart("muffleWarning")
                                } else {
                                    message(w$message)
                                }
                            })
    r <- mean(vapply(s, function(x){
        x$r.squared
    }, numeric(1)))

    r

}

