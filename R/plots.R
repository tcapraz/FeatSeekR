#' @title plotVarianceExplained
#' @description plot variance explained from 1 to max_features in res
#'
#' @param res result SummarizedExperiment from FeatSeek function
#'
#' @return returns plot of variance explained vs number of features
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(100*30), dim=c(30,100),
#'             dimnames = list(paste("feature", seq_len(30)), NULL))
#' reps <- rep(c(1,2), each=50)
#' res <- FeatSeek(data, reps, max_features=20)
#'
#' # res stores the 20 selected features ranked by their replicate reproducibility
#' plotVarianceExplained(res)
#'
#' @export
plotVarianceExplained <- function(res){
    plot(seq_len(nrow(res)), SummarizedExperiment::rowData(res)$explained_variance,
        xlab = "Number of selected features",
        ylab = "Fraction of explained variance")
}



#' @title plotSelectedFeatures
#'
#' @description plot correlation matrix of selected feature sets
#'
#' @param res result SummarizedExperiment from FeatSeek function
#' @param n_features top n_features to plot. if NULL then the maximum number
#'  of features in res will be plotted
#'
#' @return returns heatmap of selected features
#'
#' @examples
#' # run FeatSeek to select the top 20 features
#' data <-  array(rnorm(100*30), dim=c(30,100),
#'             dimnames = list(paste("feature", seq_len(30)), NULL))
#' reps <- rep(c(1,2), each=50)
#' res <- FeatSeek(data, reps, max_features=20)
#'
#' # res stores the 20 selected features ranked by their replicate reproducibility
#' # plot the top 5 features
#' plotSelectedFeatures(res, n_features=5)
#'
#' @export
plotSelectedFeatures <- function(res, n_features=NULL){
    if (is.null(n_features)){
        n_features <- nrow(res)
    }
    cor <- stats::cor(t(SummarizedExperiment::assays(res)$selected)[,seq_len(n_features)], use = "pairwise.complete.obs")

    range <- max(abs(cor))
    pheatmap::pheatmap(cor, treeheight_row = 0, treeheight_col = 0, legend=TRUE,
                    show_colnames =FALSE, show_rownames = TRUE,
                    breaks = seq(-range, range, length.out = 100))

}
