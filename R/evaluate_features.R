#' @title get_opt
#'
#' @description Get the optimal number of features based on the KS statistic
#' between the inital and the current distribution
#' of pairwise replicate correlations.
#'
#' @param res result dataframe of the FeatSeek function,
#' containing the KS statistic for each selected feature
#'
#' @return optimal number of features
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' res <- FeatSeek(data, max_features=20)
#'
#' # res stores the names, replicate correlations and
#' # the KS statistic for each selected feature
#' # we can make use of the KS statistic it to determine
#' # a smaller set of features with high degree of uniqueness
#' opt <- get_opt(res)
#' res_opt <- res[seq_len(opt),]
#'
#'
#' @export
get_opt <- function(res){
    stat <- res$ks_stat
    if (is.null(stat)) stop("Please provide output dataframe of FeatSeek with column ks_stat input!")
    if (length(stat) < 4){
        message("Less than 4 features selected! Cannot determine optimal set, returning number of selected features!")
        stop <- length(stat)
    } else{
        for (i in seq(3,length(stat))){
            # get change in ks statistic as sliding window
            window <- ifelse(i <= 6, i-1, i-5)
            if (abs(stat[[i]]-stat[[window]]) < 0.005)
                break
        }
        stop <- ifelse(i <= 6, i-2, i-6)
    }
    as.integer(stop)
}

#' @title svd_entropy
#'
#' @description Calculate SVD entropy for each selected feature set 1 ... k.
#'
#'
#' @param data input data
#' (3 dimensional array with samples x features x replicates)
#' @param res result dataframe returned by FeatSeek()
#'
#' @return list of SVD entropies for each selected feature set 1 ... k
#'
#' @export
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' res <- FeatSeek(data, max_features=20)
#'
#' # evaluate SVD entropy for selected feature set
#' svd_entropy(data, res)
svd_entropy <- function(data, res){
    # calc mean over replicates
    data.mean <- apply(data, c(1,2), mean)
    S <- res$selected
    entropy <- list()
    for (i in seq_len(length(S))){
        sel <- S[seq_len(i)]
        Sdata <- data.mean[,sel, drop=FALSE]
        singv <- svd(Sdata)[[1]]
        singv.norm <- (singv ^ 2) / (sum(singv ^ 2))
        entropy[[i]] <- -1/log(length(singv.norm )) * sum(singv.norm *log(singv.norm ))
    }
    entropy
}


#' @title reconstruction_error
#'
#' @description Calculate reconstruction error
#' for each selected feature set 1 ... k.
#'
#' @param data input data
#' (3 dimensional array with samples x features x replicates)
#' @param res result dataframe returned by FeatSeek()
#'
#' @return list of reconstruction errors for each selected feature set 1 ... k
#'
#' @export
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' res <- FeatSeek(data, max_features=20)
#'
#' # evaluate reconstruction error for selected feature set
#' reconstruction_error(data, res)
reconstruction_error <- function(data, res){
    # calc reconstruction for each selected feature set 1 ... k
    # calc mean over replicates
    data.mean <- apply(data, c(1,2), mean)
    S <- res$selected
    error <- list()
    for (i in seq_len(length(S))){
        sel <- S[seq_len(i)]
        # First we calculate the fraction of variance explained
        # by the selected feature set for each remaining feature.
        rest <- data.mean[, !dimnames(data.mean)[[2]] %in% sel, drop=FALSE]
        Sdata <- data.mean[,sel, drop=FALSE]
        # model remaining features by selected features
        model <- stats::lm(rest ~  Sdata + 0)
        # model holds a fit for each remaining feature
        # get residuals
        e <- model$residuals
        # calculate frobenius norm
        error[[i]] <- sqrt(norm(e, type="F"))
    }
    error
}
