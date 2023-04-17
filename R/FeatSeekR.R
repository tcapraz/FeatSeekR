

#' FeatSeekR an R package for unsupervised feature selection
#' 
#' FeatSeekR performs unsupervised feature selection using replicated 
#' measurements. It iteratively selects features with the highest 
#' reproducibility across conditions, after projecting out those dimensions 
#' from the data that are spanned by the previously selected features. The 
#' selected a set of features has a high replicate reproducibility and a high 
#' degree of uniqueness.
#' 
#' For information on how to use this package please type 
#' \code{vignette("FeatSeekR-vignette")}.
#' 
#' Please post questions regarding the package to the Bioconductor Support Site:
#' 
#' \url{https://support.bioconductor.org}
#' 
#' @author Tümay Capraz
#' 
#' @docType package
#' @name FeatSeekR
#' @keywords package
NULL





#' @title FeatSeek
#' @description This function ranks features of a 2
#' dimensional array according to their reproducibility between conditions.
#'
#'
#' 
#' @param data \code{SummarizedExperiment} with assay named \code{data}, 
#' where samples
#' belongs to different conditions. Which sample belongs to which condition
#' should be indicated in \code{colData} slot conditions. 
#' Or \code{matrix} with features x samples.
#' Each conditions have multiple samples from replicated measurements.
#'
#' @param conditions factor of length samples,
#' indicating which sample belongs to which condition. Only required if 
#' \code{data} is provided as \code{matrix}.
#' @param init \code{character vector} with names of initial features.
#' If \code{NULL} the feature with highest F-statistic will be used
#' @param max_features \code{integer} number of features to rank
#' @param verbose \code{logical} indicating whether messages should be printed
#' 
#' @return \code{SummarizedExperiment} containing one assay with the 
#' selected features. \code{rowData} stores for each selected feature the 
#' F-statistic under \code{metric},
#' the cumulative explained variance under \code{explained_variance} and
#' the feature names under \code{selected}
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(100*30), dim=c(30, 100),
#' dimnames <- list(paste("feature", seq_len(30)), NULL))
#' conds <- rep(seq_len(50), 2)
#' res <- FeatSeek(data, conds, max_features=20)
#'
#' # res stores the 20 selected features ranked by their replicate 
#' # reproducibility
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @export
FeatSeek <- function(
    data, 
    conditions=NULL, 
    max_features=NULL, 
    init=NULL, 
    verbose=TRUE) {
    stopifnot(
        is.numeric(conditions) | is.null(conditions),
        is.numeric(max_features) | is.null(max_features),
        is.character(init) | is.null(init),
        is.logical(verbose)
    )
    se <- check_input(data, max_features, conditions)

    p <- nrow(se)
    n <- ncol(se)
    c <- length(unique(se$conditions))
    r <- n/c

    # initialize starting set of features
    init <- init_selected(init, se)
    if (verbose){
        message("Input data has: \n",
                n, " samples \n",
                c, " conditions \n",
                p, " features \n",
                r, " replicates")
    }
    # init max_features
    if(is.null(max_features)){
        max_features <- p
    }
    # initialize output dataframe
    res <-
        data.frame(
            metric = rep(NA, max_features),    # selection metric
            selected = rep("", max_features),  # name of selected features
            explained_variance = rep(NA, max_features)
            # variance explained by selected features
        )
    metric_all <- list()
    
    # initialize matrix of selected features
    S <- array(NA, dim=c(n, p))
    k <- 1
    data <- t(assay(se, "data"))
    conditions <- se$conditions
    data0 <- data
    # start ranking features
    # in each iteration we select the feature with
    # highest consistency between conditions
    # after projecting out the previously selected ones
    if (verbose){
        message("Starting feature ranking!")
    }
    while (k <= max_features)  {
        if (k > 1) {
            data <- fit_lm(data, S, k)
        }
        # check if any replicate contains only zero residuals
        if (any(apply(data, 2, function(x) {all(x == 0)}))) {
            # if all residuals are close to zero stop selection procedure
            if (verbose){
                message(
                    "Stopping selection procedure with k=",
                    k - 1,
                    " features because remaining features are linear 
                    combinations of already selected features!"
                )
            }
            break
        }
        # calculate mean pairwise correlations between all conditions
        metric <- calcFstat(data, conditions)
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
        S[, k] <- data[, I, drop = FALSE]
        # drop selected feature from data
        data <- data[, dimnames(data)[[2]]!= I, drop = FALSE]
        if (verbose){
            message("Iteration: ", k, " selected = ", res$selected[k],
                    ", replicate F-statistic = ", res$metric[k])
        }
        k <- k + 1
    }
    if (verbose){
        message("Finished feature ranking procedure!")
        message("Calculating explained variance of selected feature sets!")
    }
    res <- res[seq_len(k-1), ]
    res$explained_variance <- vapply(seq_along(res$selected), function(i){
        variance_explained(data0, res$selected[seq_len(i)])
    }, numeric(1))
    se_res <- SummarizedExperiment(
                assays=list(selected=t(data0[, res$selected, drop=FALSE])), 
                rowData=res
                )
    se_res
}




#' @title calcFstat
#'
#' @param data 2 dimensional \code{array} with samples x features, 
#' where samples
#' belongs different conditions. The function was adapted from the
#' \code{genefilter} package.
#' @param fac \code{factor} of length samples,
#' indicating which sample belongs to which condition
#' 
#' @return F-statistic for all features
#'
#' @importFrom methods is
#' 
#' @keywords internal
calcFstat <- function(data, fac){
    stopifnot(
        length(fac)==nrow(data), 
        is.factor(fac), 
        is.matrix(data)
    )
    data <- scale(data)
    data <- t(data)
    data   <- data[, !is.na(fac), drop=FALSE]
    fac <- fac[!is.na(fac)]
    
    # Number of levels (groups)
    k <- nlevels(fac)
    
    # group_m: a ncol(data) x nlevels(fac) matrix with the means of each 
    # condition
    group_m <- matrix(
        vapply(levels(fac), function(fl) 
            rowMeans(data[, which(fac==fl), drop=FALSE]), numeric(nrow(data))),
        nrow = nrow(data),
        ncol = nlevels(fac))
    
    # group_gm: a matrix of condition means, with as many rows as data, 
    # columns correspond to conditions
    group_m <- group_m[, fac, drop=FALSE]
    
    # degree of freedom 1
    dff    <- k - 1
    
    # all_m: a matrix of same size as data with overall means
    all_m <- matrix(rowMeans(data), ncol=ncol(data), nrow=nrow(data))
    
    # degree of freedom 2
    dfr    <- ncol(data) - dff
    
    # mean sum of squares
    mssf   <- rowSums((all_m - group_m)^2) / dff
    mssr   <- rowSums((data - group_m)^2) / dfr
    
    # F statistic
    fstat  <- mssf/mssr
    fstat
}
