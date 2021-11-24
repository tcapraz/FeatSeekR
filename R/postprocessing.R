#' @title stop_criterion
#'
#' @description Get the optimal number of features based on the KS statistic
#' between the inital and the current distribution of pairwise replicate correlations.
#'
#' @param dist distribution of pairwise replicate pearson correlations
#'
#' @return returns optimal number of features
#'
#' @keywords internal
stop_criterion <- function(stat){
  for (i in seq(3,length(stat))){
    # get change in ks statistic as sliding window
    window <- ifelse(i <= 6, i-1, i-5)
    if (abs(stat[[i]]-stat[[window]]) < 0.005)
      break
  }
  stop <- ifelse(i <= 6, i-2, i-6)
  stop
}


#'
#'
#' #' @title Postprocessing
#' #' @description  After the selection procedure we determine the optimal minimum set of features.
#' #' The optimal set maximizes the fraction of variance of the remaining features
#' #' that can be explained by the selected set of features, while minimizing redundancy.
#' #' @param data 3 dimensional array with samples x replicates x features.
#' #' @param result result dataframe of feature selection run with FeatSeek
#' #' @param lambda regulates tradeoff between variance explained and redundancy
#' #'
#' #' @return list with original dataframe with variance explained and redundancy metrics for each selected feature added and optimal number of features
#' #' @export
#' #'
#' postprocessing <- function(data, result, lambda=0.5){
#'
#'     n <- dim(data)[1]
#'     r <- dim(data)[2]
#'     p <- dim(data)[3]
#'
#'     # check if data has feature names
#'     if(!is.null(dimnames(data)[[3]])) {
#'         features <- dimnames(data)[[3]]
#'     }else {
#'         stop("No feature names given or features not in correct dimension of data array!")
#'     }
#'
#'     if (!"selected" %in% names(result)){
#'         stop("No selected features in selection  result!\n Please pass output of FeatSeek as result argument.")
#'     }
#'     if (length(result$selected) < 2 ){
#'         stop("Number of selected features must be >= 2!")
#'     }
#'
#'     if (any(!result$selected %in% dimnames(data)[[3]])){
#'         stop("Selected features not in input data!")
#'     }
#'     # average over replicates
#'     mean_data <- apply(data, MARGIN=c(1,3), sum)/r
#'     # selected features
#'     selected_data <- mean_data[,result$selected]
#'
#'     result$var_explained <- rep(NA, dim(result)[1])
#'     result$redundancy <- rep(NA, dim(result)[1])
#'
#'     # First we calculate the fraction of variance explained by the selected feature set for each remaining feature.
#'     for (i in seq_len(length(result$selected))){
#'         sel <- selected_data[,result$selected[seq_len(i)]]
#'         rest <- mean_data[, !dimnames(mean_data)[[2]] %in% result$selected[1:i], drop=FALSE]
#'         # model remaining features by selected features
#'         model <- lm(rest ~  sel + 0)
#'         # model holds a fit for each remaining feature
#'         summary <- summary(model)
#'         # get average of adjusted r square values over all remaining features
#'         r_square <- mean(unlist(lapply(summary, function(x) {x$adj.r.squared})))
#'         result$var_explained[[i]] <- r_square
#'     }
#'
#'     # Next we check for redundancy in all possible feature sets of our feature ranking. T
#'     # We measure redundancy by SVD entropy. SVD entropy is defined as the normalized entropy of the singular values of a SVD of all selected features.
#'     # A feature set with low SVD entropy has higher redundancy as most of the variance can be captured by few singualar values.
#'
#'     entropy <- rep(NA, length(result$selected))
#'     for (i in seq_len(length(result$selected))){
#'         sel <- selected_data[,result$selected[seq_len(i)]]
#'         singv <- svd(sel)[[1]]
#'         singv.norm <- (singv^2) / (sum(singv^2))
#'         entropy[i] <- calc_entropy(singv.norm)
#'     }
#'     entropy[1] <- 1
#'     result$redundancy <- entropy
#'     # get optimum number of selected features which maximizes both variance explained of remaining features and SVD entropy
#'     # lambda controls the contribution of redundancy and variance explained
#'     opt <- which.max((result$var_explained*lambda) + result$redundancy*(1-lambda))
#'     # return optimal set
#'     list(result=result, opt_features = opt)
#' }
