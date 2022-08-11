




#' @title FeatSeek
#' @description This function ranks features of a 3
#' dimensional array according to their consistency between replicates.
#'
#'
#'
#' @param data 3 dimensional array with samples x features x replicates
#' @param init vector with names of initial features.
#' If NULL the feature with highest replicate correlation will be used
#' @param max_features integer number of features to rank
#' @param filter_thr float between 0 and 1.
#' Features with correlation < filter_thr are removed
#' @return Dataframe with selected feature names,
#' replicate correlation and KS statistic
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' res <- FeatSeek(data, max_features=20)
#'
#' # res stores the 20 selected features ranked by their replicate reproducibility
#'
#' @export
FeatSeek <- function(data,  max_features, init=NULL, filter_thr = NULL) {

    data <- check_input(data)
    n <- dim(data)[1]
    p <- dim(data)[2]
    r <- dim(data)[3]

    # initialize starting set of features
    init <- init_selected(init, data)

    message("Input data has: \n",
            n, " samples \n",
            r, " replicates \n",
            p, " features")

    # remove features with correlation < filter_thr
    if(!is.null(filter_thr)){
        data <- filter(data, init, filter_thr)
    }
    n <- dim(data)[1]
    p <- dim(data)[2]
    r <- dim(data)[3]
    # initialize output dataframe
    res <-
        data.frame(
            metric = rep(NA, max_features),    # selection metric
            selected = rep("", max_features),  # name of selected features
            explained_variance = rep(NA, max_features) # variance explained by
                                                       # selected features
        )

    metric_all <- list()

    # initialize matrix of selected features
    S <- array(NA, dim=c(n,p,r))
    k <- 1

    data0 <- data
    # start ranking features
    # in each iteration we select the feature with
    # highest consistency between replicates
    # after projecting out the previously selected ones
    message("Starting feature ranking!")
    while (k <= max_features)  {
        if (k > 1) {
            data <- fit_lm(data, S, k)
        }
        # check if any replicate contains only zero residuals
        if (any(apply(data, 3, function(x) {all(x==0)}))) {
            # if all residuals are close to zero stop selection procedure
            message(
                "Stopping selection procedure with k=",
                k - 1,
                " features because remaining features are linear combinations
                of already selected features!"
            )
            break
        }
        # calculate mean pairwise correlations between all replicates
        metric <- calcFstat(data)
        metric_all[[k]] <- metric
        names(metric) <- dimnames(data)[[2]]
        if (k > length(init)) {
            # select feature whose residuals have the highest correlation
            I <- names(metric)[which.max(metric)]
        } else{
            # first features are set to init ones
            I <- init[k]
        }
        # store metric of selected feature
        res$metric[k] <- metric[I]
        # store selected feature
        res$selected[k] <- I
        S[, k,] <- data[, I, , drop = FALSE]
        # drop selected feature from data
        data <- data[, dimnames(data)[[2]]!= I,  , drop = FALSE]
        message("Iteration: ",k," selected = ",res$selected[k],
                ", replicate F-statistic = ",res$metric[k])
        k <- k + 1

    }
    message("Finished feature ranking procedure!")
    res <- res[seq_len(k-1),]

    message("Calculating explained variance of selected feature sets!")
    res$explained_variance <- vapply(seq_along(res$selected), function(i){
        variance_explained(data0,res$selected[seq_len(i)])
    }, numeric(1))
    res
}

#' @title calcCor
#'
#' @description Calculate pairwise replicate correlations
#'
#' @param data 3 dimensional array with samples x features x replicates
#' @param r number of replicates
#'
#' @import matrixStats
#' @import rstatix
#'
#'
#' @return pairwise replicate pearson correlations
#'
#'
#' @keywords internal
calcCor <- function(data, r){
    feat_cor <- stats::cor(data, use="pairwise.complete.obs")
    feat_cor[is.na(feat_cor)] <- 0
    feat_cor[upper.tri(!diag(nrow=r))]
}


#' calcFstat
#'
#' @description Calculate F statistic between replicates
#'
#' @param data 3 dimensional array with samples x features x replicates
#' @param scale logical indicating whether to center and scale the data
#' @param complete logical indicating whether the data is complete or samples
#'                are missing
#'
#' @return  F-like statistic
#'
#'@keywords internal
calcFstat <- function(data, scale=TRUE, complete=TRUE){
    # scale each replicate
    if (scale){
        data <- lapply(seq_len(dim(data)[3]),  function(x, data){
            scale(data[,,x],center=TRUE, scale=TRUE)
            }, data=data)
        data <- abind::abind(data, along=3)

    }
    Fvalue <- vapply(seq_len(dim(data)[2]), function(x, data, complete){
        if (complete) {
            v <- data[,x,]
            varwithin  <- mean(matrixStats::rowVars(t(v)))
            varbetween <- var(rowMeans(v))
            sqrt(varbetween / varwithin)
        } else {
            v <- data[,x,]
            grpsb <- c(vapply(seq_len(dim(v)[2]), function(x, v){
                seq_along(v[,x])
            }, integer(nrow(v)), v=v))
            grpsw <- c(vapply(seq_len(dim(v)[1]), function(x, v){
                seq_along(v[x,])
            }, integer(ncol(v)), v=v))
            vt <- t(v)
            varwithin  <- tapply(vt, grpsw, FUN = stats::var,  na.rm = TRUE) |> mean()
            varbetween <- tapply(v, grpsb, FUN = mean, na.rm = TRUE) |> stats::var()
            sqrt(varbetween / varwithin)
        }

    },FUN.VALUE = numeric(1),  data=data, complete=complete)
    Fvalue
}

