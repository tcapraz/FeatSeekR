#' @title calc_metric
#'
#' @description Calculate pairwise replicate correlations
#'
#' @param data 3 dimensional array with samples x features x replicates
#' @param r number of replicates
#'
#' @return pairwise replicate pearson correlations
#'
#'
#' @keywords internal
calc_metric <- function(data, r){
    feat_cor <- stats::cor(data, use="pairwise.complete.obs")
    feat_cor[is.na(feat_cor)] <- 0
    feat_cor[upper.tri(!diag(nrow=r))]
}


#' @title FeatSeek
#' @description This function ranks features of a 3 dimensional array according to their consistency between replicates.
#'
#'
#'
#' @param data 3 dimensional array with samples x features x replicates.
#' @param init vector with names of initial features. If NULL the feature with highest replicate correlation will be used.
#' @param max_features integer number of features to rank
#' @param filter_thr float between 0 and 1. Features with correlation < filter_thr are removed .
#' @return Dataframe with selected features, rmsd, 1/VIF and variance explained for the selection process
#'
#' @export
FeatSeek <- function(data,  max_features, init=NULL, filter_thr = NULL) {

    data <- check_input(data)

    n <- dim(data)[1]
    p <- dim(data)[2]
    r <- dim(data)[3]

    # initialize starting set of features
    init <- init_selected(init, data)

    message("Input data has: \n", n, " samples \n", r, " replicates \n", p, " features")

    # remove features with correlation < filter_thr
    if(!is.null(filter_thr)){
        data <- filter(data, init, filter_thr)
    }

    # initialize output dataframe
    res <-
        data.frame(
            metric = rep(NA, max_features),    # selection metric for each iteration
            selected = rep("", max_features),  # name of selected features
            ks_stat = rep(NA, max_features)
        )

    # calculate initial distribution of correlations for all features and replicates
    initial_dist <- as.vector(apply(data, 2, calc_metric, r=r))

    # initialize matrix of selected features
    S <- matrix(NA, nrow=n, ncol=max_features)
    k <- 1
    nonzero_residuals <- TRUE
    continue <- TRUE

    # start ranking features
    # in each iteration we select the feature with highest consistency between replicates
    # after projecting out the previously selected ones
    message("Starting feature ranking!")
    while (k <= max_features & continue==TRUE & nonzero_residuals ==TRUE)  {
        if (k > 1) {
            data <- fit_lm(data, S, k, r)
        }
        if(all(data==0)) {
            # if all residuals are close to zero stop selection procedure
            message(
              "Stopping selection procedure with k=",
              k - 1,
              " features because remaining features are linear combinations of already selected features. Consider lowering max_features!"
            )
            break
        }
        # calculate mean pairwise correlations between all replicates
        allmetrics <- apply(data, 2, calc_metric, r=r)
        metric <-   colMeans(allmetrics)

        # KS statistic to compare distribution of metrics to previous iteration
        ks <- stats::ks.test(as.vector(allmetrics), initial_dist)
        res$ks_stat[k] <- ks$statistic

        if (k > length(init)) {
            names(metric) <- dimnames(data)[[2]]
            # select feature whose residuals have the highest correlation
            I = names(metric)[which.max(metric)]
        } else{
            # first features are set to init ones
            I = init[k]
        }
        # store metric of selected feature
        res$metric[k] = metric[I]
        # store selected feature
        res$selected[k] <- I
        # add mean between replicates to selected subset
        S[, k] = apply(data[, I, , drop = FALSE], 1, mean, na.rm = TRUE)

        # drop selected feature from data
        data = data[, dimnames(data)[[2]]!= I,  , drop = FALSE]
        message("Iteration: ",k," selected = ",res$selected[k],", replicate correlation = ",res$metric[k], ", KS statistic = ", round(res$ks_stat[k], digits=2))
        k <- k + 1
    }
    message("Finished feature ranking procedure!")
    res <- res[seq_len(k-1),]
}
