#' #' @title reconstruction_error
#' #'
#' #' @description Calculate reconstruction error
#' #' for each selected feature set 1 ... k.
#' #'
#' #' @param data input data
#' #' (3 dimensional array with samples x features x replicates)
#' #' @param res result dataframe returned by FeatSeek()
#' #'
#' #' @return list of reconstruction errors for each selected feature set
#' 1 ... k
#' #'
#' #' @export
#' #' @examples
#' #' # run FeatSeek to select the top 20 features
#' #' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' #' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' #' res <- FeatSeek(data, max_features=20)
#' #'
#' #' # evaluate reconstruction error for selected feature set
#' #' reconstruction_error(data, res)
#' reconstruction_error <- function(data, res){
#'     # calc reconstruction for each selected feature set 1 ... k
#'     # calc mean over replicates
#'     data.mean <- apply(data, c(1,2), mean)
#'     y <- data
#'     S <- res$selected
#'     X <- data[,S]
#'
#'     P = X %*% inv(t(X)%*%X) %*%t(X)
#'     y_ <- apply(y,2, function(x)x %*% P)
#'     e <- y-y_
#'
#'     # calculate frobenius norm
#'     if (length(dim(e)< 2)){
#'       error <- norm(e, type="F")
#'     } else{
#'       error <- norm(e, type="2")
#'     }
#'     sqrt(error)
#' }
#'



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

