
calc_metric <- function(data,r){
  feat_cor <- cor(data, use="pairwise.complete.obs")
  feat_cor[is.na(feat_cor)] <- 0
  negative <-  any(feat_cor < 0)
  res <-   mean(feat_cor[upper.tri(!diag(nrow=r))])
  if (negative == TRUE){
    res <- sign(res)*-1*res
  }
  res

}

calc_metric_all <- function(data,r){
  feat_cor <- cor(data, use="pairwise.complete.obs")
  feat_cor[is.na(feat_cor)] <- 0

  res <- feat_cor[upper.tri(!diag(nrow=r))]

  res
}

filter <- function(data, filter_thr){
  message("Filtering out features with low initial correlation across replicates....")
  # drop features with sd == 0
  sdl <- lapply(seq_len(dim(data)[2]), function(i, data) {
    apply(data[, i, ], 2, function(x)
      sd(x))
  }, data = data)
  zeros <- lapply(sdl, function(i){
    i == 0
  })
  zeros <- Reduce("|", zeros)
  data <- data[,,!zeros]
  # get mean of pairwise correlation between replicates, i.e. the off-diagonal values of the correlation matrix
  cor <- apply(data,3, function(x) mean(cor(x, use="pairwise.complete.obs")[upper.tri(!diag(nrow=dim(data)[2]))]))
  keep <-  cor_raw > filter_thr
  features <- dimnames(data)[[3]]
  # check if preselected features would be filtered out and issue a warning
  if(any(!preselected %in% features[keep]))
    warning("The following preselected features have replicate correlation lower than filter_thr:\n",
            paste(preselected[!preselected %in% features[keep]],collapse=","),
            "\nConsider removing them!")
  features <- union(features[keep], preselected)
  features <- features[!is.na(features)]
  data <- data[, , features]
  p <- dim(data)[3] # get number of remaining features
  message(paste("Remaining features:", p))
  data
}


#' @title FeatSeek algorithm
#' @name FeatSeek
#' @description This function ranks features of a 3 dimensional array according to their consistency between replicates.
#'
#'
#'
#' @param data 3 dimensional array with samples $\times$ replicates x features.
#' @param preselected vector with names of preselected features.
#' @param max_features integer number of features to rank
#' @param filter_thr float between 0 and 1. Features with correlation < filter_thr are removed .
#' @return Dataframe with selected features, rmsd, 1/VIF and variance explained for the selection process
#'
#' @export
FeatSeek <- function(data,  max_features, preselected=NULL, filter_thr = NULL, stop_thr=0.5, ret_data=FALSE) {

  n <- dim(data)[1]
  r <- dim(data)[2]
  p <- dim(data)[3]

  message(paste0("Starting feature selection!"))
  message(paste0("Input data has: \n", n, " samples \n", r, " replicates \n", p, " features"))

  # check and filter inputs
  if(!r>=2) stop("At least 2 replicates required!")

  # check if data has feature names
  if(!is.null(dimnames(data)[[3]])) {
    features <- dimnames(data)[[3]]
  }else {
    stop("No feature names given or features not in correct dimension of data array!")
  }

  if(!is.null(preselected) & !all(preselected %in% features)){
    stop("Could not find preselected features in data!")
  }
  if (is.null(preselected)){
    cor <- apply(data,3, function(x) mean(cor(x, use="pairwise.complete.obs")[upper.tri(!diag(nrow=dim(data)[2]))]))
    names(which.max(cor))
  }
  # remove features with correlation < filter_thr
  if(!is.null(filter_thr)){
    data <- filter(data, filter_thr)
  }

  # initialize output dataframe
  res <-
    data.frame(
      metric = rep(NA, max_features),    # selection metric for each iteration
      selected = rep("", max_features),  # name of selected features
      ks_stat = rep(NA, max_features),
      rec_error = rep(NA, max_features),
      svd_entropy = rep(NA, max_features)
    )
  allmetrics <- list()
  alldata <- list()
  # initialize matrix of selected features
  sel <- matrix(NA, nrow=n, ncol=max_features)
  k <- 1
  nonzero_residuals <- TRUE
  continue <- TRUE

  # start ranking features
  # in each iteration we select the feature with highest consistency between replicates after projecting out the previously selected ones
  while (k <= max_features & continue==TRUE & nonzero_residuals ==TRUE)  {
    if (k > 1) {
      # get current features in data
      features <- dimnames(data)[[3]]
      d <- list()
      # fit linear model for each replicate and overwrite data with residuals
      s <- sel[,1:(k-1), drop = FALSE]
      for (j in seq_len(r)){
        #s <- sel[[r]][,1:(k-1), drop = FALSE]
        model = lm(data[, j, ] ~ s + 0, na.action = na.exclude)
        resids <- array(NA, dim(data[, j,]))
        resids[!is.na(data[, j,])] <- model$residuals
        data[, j, ] = ifelse(abs(resids)<10^(-10), 0, resids)
      }
    }
    alldata[[k]] <- data
    if(!all(data==0)) {

      # calculate mean pairwise correlations between all replicates
      # set to negative if one pair is negatively correlated
      metric <- apply(data, 3, calc_metric, r=r)
      allmetrics[[k]] <- as.vector(apply(data, 3, calc_metric_all, r=r))
      # KS statistic to compare distribution of metrics to previous iteration
      if (k>1){
        metric_std <- scale(metric)

        ks <- ks.test(metric_std, pnorm)

        res$ks_stat[k] <- ks$statistic
        res$ks_pval[k] <- ks$p.value

      }


      if (k > length(preselected)) {
        names(metric) <- dimnames(data)[[3]]
        # select feature whose residuals have the highest correlation
        I = names(metric)[which.max(metric)]
      } else{
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
      if (ratio_positive <= stop_thr) {
        continue <- FALSE
        message(paste0("Stopping selection procedure as ratio of positive correlations <= 0.5!"))
      }
      k <- k + 1
    }else {
      # if all residuals are close to zero stop selection procedure
      message(paste0("Stopping selection procedure with k=",k-1," features because remaining features are linear combinations of already selected features. Consider lowering max_features!"))
      nonzero_residuals <- FALSE
    }
  }
  if (ret_data ==TRUE){
    return(list(res[1:(k-1),], allmetrics, alldata, sel))
  }
  else{
    return(list(res[1:(k-1),], allmetrics))
  }
}
