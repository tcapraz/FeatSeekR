


#' @title FeatSeek algorithm
#' @name FeatSeek
#' @description This function ranks features of a 3 dimensional array according to their consistency between replicates.
#' After the ranking procedure the optimal subset of features is determined by maximizing the inverse of the variance inflation factor (VIF)
#' and the fraction of variance of the remaining features that can be explained by the selected ones. The VIF is a measure for collinearity
#' in linear regression, the higher the inverse of VIF the lower the redundancy.
#'
#' @param data 3 dimensional array with samples $\times$ replicates x features.
#' @param preselected vector with names of preselected features.
#' @param max_p integer number of features to rank
#' @param lambda  float between 0 and 1. Setting the tradeoff between variance explained and redundancy
#'                if lambda = 1 only the consistency of residuals between replicates is considered feature selection.
#' @param filter_thr float between 0 and 1. Features with correlation < filter_thr are removed .
#' @param threads int that specifies how many cores to use
#'
#' @return Dataframe with selected features, rmsd, 1/VIF and variance explained for the selection process
#'
#' @export
FeatSeek <- function(data, preselected, max_p, filter_thr = 0.5, lambda = 0.5) {

  n <- dim(data)[1]
  r <- dim(data)[2]

  # check and filter inputs
  if(!r>=2) stop("At least 2 replicates required!")

  # check if data has feature names
  if(!is.null(dimnames(data)[[3]])) {
    features <- dimnames(data)[[3]]
  }else {
    stop("No feature names given or features not in correct dimension of data array!")
  }

  if(!all(preselected %in% features)){
    stop("Could not find preselected features in data!")
  }

  # remove features with correlation < filter_thr
  if(filter_thr>0){
    message("Filtering out features with low initial correlation across replicates....")
    # get mean of pairwise correlation between replicates, i.e. the off-diagonal values of the correlation matrix
    cor_raw <- apply(data,3, function(x) mean(cor(x, use="pairwise.complete.obs")[!diag(TRUE,r)]))
    include <-  which(cor_raw > filter_thr)
    # check if preselected features would be filtered out and issue a warning
    if(any(!preselected %in% features[include]))
      warning("The following preselected features have replicate correlation lower than filter_thr:\n",
              paste(preselected[preselected %in% features[include]],collapse=","),
              "\nConsider removing them!")
    features <- union(features[include], preselected)
    data <- data[, , features]
    p <- dim(data)[3] # get number of features
    remaining_feat <- dimnames(data)[3]
    message(paste("Remaining features:", p))
  }
  # remove correlation matrix from memory
  rm(cor_raw)

  # store mean over all replicates of original data
  mean_data <- apply(data, MARGIN=c(1,3), sum)/r

  # initialize output dataframe
  res <-
    data.frame(
      metric = rep(NA, max_p),    # selection metric for each iteration
      selected = rep("", max_p),  # name of selected features
      vif = rep(NA, max_p),       # inverse of variance inflation factor
      r_square = rep(NA, max_p),  # variance explained of remaining features
      ratio_positive = rep(NA, max_p)
    )
  # initialize matrix of selected features
  sel <- matrix(NA, nrow=n, ncol=max_p)
  k <- 1
  nonzero_residuals <- TRUE
  continue <- TRUE

  # start ranking features
  # in each iteration we select the feature with highest consistency between replicates after projecting out the previously selected ones
  while (k <= max_p & continue==TRUE & nonzero_residuals ==TRUE)  {
    if (k > 1) {
      # get current features in data
      features <- dimnames(data)[[3]]
      d <- list()
      # fit linear model for each replicate and overwrite data with residuals
      # the data is split into chunks of size chunksize and each thread fits
      for (j in seq_len(r)){
        s <- sel[,1:(k-1), drop = FALSE]
        d[[j]] <-  abind(applyKernel(data[,j,], fit_lm, threads, sel=s ), along=3)
      }
      # concatenate to 3 dimensional array
      data <- abind(d, along=2)
      # set feature names, as the dimnames are lost when fitting the lm
      dimnames(data)[[3]] <- features
    }
    if(!all(data==0)) {
      # calculate rmsd of residuals between replicates, as a measure of consistency
      # metric = apply(data, 3, function(x) {
      #   # normalize each residual vector
      #   sums <- sqrt(colSums(x ^ 2))
      #   X <- sweep(x, 2, sums, "/")
      #   # get mean vector between replicates
      #   m = rowMeans(X)
      #   # calculate pairwise rmsd to the mean vector
      #   rmsd = apply(X,2,function(x2,m){
      #     sqrt(mean((m-x2)^2))
      #   }, m=m)
      #   # sum over the replicates
      #   sum(rmsd)
      # })
      metric <- apply(data, 3, function(x){
        feat_cor <- cor(x, use="pairwise.complete.obs")
        negative <-  any(feat_cor < 0)
        res <-   mean(feat_cor[!diag(TRUE,r)])
        if (negative == TRUE){
          res <- res*-1
        }
        res
      })

      if (k > length(preselected)) {
        names(metric) <- dimnames(data)[[3]]
        # select feature whose residuals have the smallest rmsd
        I = names(metric)[which.max(metric)]
      }
      if (k <= length(preselected)) {
        # first features are set to preselected ones
        I = preselected[k]
      }
      # store metric of selected feature
      res$metric[k] = metric[I]
      # store selected feature
      res$selected[k] <- I
      # add mean between replicates to selected subset
      sel[, k] = apply(data[, , I, drop = FALSE], 1, mean, na.rm = TRUE)
      # drop selected feature from data
      data = data[, , dimnames(data)[[3]] != I, drop = FALSE]
      ratio_positive <- sum(metric>0)/length(metric)
      res$ratio_positive[k] <- ratio_positive
      message("k=",k," selected = ",res$selected[k],", metric = ",res$metric[k], ", ratio positive = ", res$ratio_positive[k])
      ratio_positive <- sum(metric>0)/length(metric)
      if (ratio_positive <= 0.5) {
        continue <- FALSE
        message(paste0("Stopping selection procedure as ratio of positive correlations <= 0.5!"))
      }
      k <- k + 1
    }else {
      # if all residuals are close to zero stop selection procedure
      message(paste0("Stopping selection procedure with k=",k-1," features because remaining features are linear combinations of already selected features. Consider lowering max_p!"))
      nonzero_residuals <- FALSE
    }
  }

  return(res)
}
