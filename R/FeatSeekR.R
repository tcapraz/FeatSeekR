#' @title FeatSeek
#' @description This function ranks features of a 2
#' dimensional array according to their reproducibility between replicates.
#'
#'
#'
#' @param data 2 dimensional array with samples x features, where each sample
#' belongs to a different replicate
#'
#' @param replicates numeric vector of length samples,
#' indicating which sample belongs to which replicate
#' @param init vector with names of initial features.
#' If NULL the feature with highest F-statistic will be used
#' @param max_features integer number of features to rank
#' @return Dataframe with selected feature names,
#' F-statistic and explained variance
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(100*30), dim=c(100,30),
#' dimnames <- list(NULL, paste("feature", seq_len(30))))
#' reps <- rep(c(1,2), each=50)
#' res <- FeatSeek(data, reps, max_features=20)
#'
#' # res stores the 20 selected features ranked by their replicate reproducibility
#'
#' @export
FeatSeek <- function(data, replicates, max_features=NULL, init=NULL) {

    check_input(data, replicates, max_features)
    n <- dim(data)[1]
    p <- dim(data)[2]
    r <- length(unique(replicates))

    # initialize starting set of features
    init <- init_selected(init, data, replicates)

    message("Input data has: \n",
            n, " samples \n",
            r, " replicates \n",
            p, " features")

    n <- dim(data)[1]
    p <- dim(data)[2]
    r <- length(unique(replicates))

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
    S <- array(NA, dim=c(n,p))
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
        if (any(apply(data, 2, function(x) {all(x==0)}))) {
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
        metric <- calcFstat(data, replicates)
        metric_all[[k]] <- metric
        names(metric) <- dimnames(data)[[2]]
        if (k > length(init)) {
            # select feature whose residuals have the highest correlation
            I <- names(metric)[which.min(metric)]
        } else{
            # first features are set to init ones
            I <- init[k]
        }
        # store metric of selected feature
        res$metric[k] <- metric[I]
        # store selected feature
        res$selected[k] <- I
        S[, k] <- data[, I,  drop = FALSE]
        # drop selected feature from data
        data <- data[, dimnames(data)[[2]]!= I,   drop = FALSE]
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



#' @title calcFstat
#'
#' @param data 2 dimensional array with samples x features, where each sample
#' belongs to a different replicate
#' @param reps numeric vector of length samples,
#' indicating which sample belongs to which replicate
#' @param scale whether to scale the data, default = TRUE
#'
#' @return F-statistic for all features
#'
#' @internal
#'
calcFstat <- function(data,reps, scale=TRUE){

    data <- scale(data)
    f <- vapply(seq_len(dim(data)[2]), function(x){
        m <- stats::lm(data[,x] ~as.factor(reps))
        s <- summary(m)
        f <- s$fstatistic[1]
        }, numeric(1))
    f
}
