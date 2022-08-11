#' @title plotVarianceExplained
#' @description plot variance explained from 1 to max_features in res
#'
#' @param res result dataframe from FeatSeek function
#'
#' @return returns plot of variance explained vs number of features
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' res <- FeatSeek(data, max_features=20)
#'
#' # res stores the 20 selected features ranked by their replicate reproducibility
#' plotVarianceExplained(res)
#'
#' @export
plotVarianceExplained <- function(res){
    plot(seq_len(nrow(res)), res$explained_variance,
        xlab = "Number of selected features",
        ylab = "Fraction of explained variance")
}



#' @title plotSelectedFeatures
#'
#' @description plot correlation matrix of selected feature sets
#'
#' @param data 3 dimensional array with samples x features x replicates
#' @param res result dataframe from FeatSeek function
#' @param n_features number of features to plot. if NULL then the maximum number
#'  of features in res will be plotted
#'
#' @return returns heatmap of selected features
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(50*30*2), dim=c(50,30,2),
#' dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
#' res <- FeatSeek(data, max_features=20)
#' # res stores the 20 selected features ranked by their replicate reproducibility
#' # plot the top 5 features
#' plotSelectedFeatures(data, res, n_features=5)
#'
#' @export
plotSelectedFeatures <- function(data, res, n_features=NULL){
  if (length(dim(data)) == 3){
    cnames <- dimnames(data)[[2]]
    data <- apply(data, c(1,2), mean )
    colnames(data) <- cnames
  }
  if (is.null(n_features)){
    n_features <- nrow(res)
  }
  cor <- stats::cor(data[,res$selected[seq_len(n_features)]], use = "pairwise.complete.obs")

  range <- max(abs(cor))
  p <- pheatmap::pheatmap(cor, treeheight_row = 0, treeheight_col = 0, legend=TRUE,
                     show_colnames =FALSE, show_rownames = TRUE,
                     breaks = seq(-range, range, length.out = 100))
  p
}
