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
#' @importFrom SummarizedExperiment rowData
#' @importFrom methods is
#' @export
plotVarianceExplained <- function(res){
    stopifnot(
        is(res, "SummarizedExperiment")
    )
    if (is.null(rowData(res)$explained_variance)){
        stop("'res' does not contain attribute explained_variance 
        Please run FeatSeek first!")
    }

    plot(seq_len(nrow(res)), res$explained_variance,
        xlab="Number of selected features",
        ylab="Fraction of explained variance")
}



#' @title plotSelectedFeatures
#'
#' @description plot correlation matrix of selected feature sets
#'
#' @param res result SummarizedExperiment from FeatSeek function
#' @param n_features top n_features to plot. if NULL then the maximum number
#'  of features in res will be plotted
#' @param assay assay slot to plot from result SummarizedExperiment object, 
#'  default is the selected features slot
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
#' @importFrom stats cor
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom pheatmap pheatmap
#' @importFrom methods is
#' 
#' @export
plotSelectedFeatures <- function(res, n_features=NULL, assay="selected"){
    stopifnot(
        is(res, "SummarizedExperiment"),
        # 'n_features' should be a scalar integer and 
        # mustn't exeed the number of available features
        is.numeric(n_features) | is.null(n_features), 
        # assay should be character that matches a unique assay name
        is.character(assay), length(assay) == 1,
        sum(grepl(assay, assayNames(res))) == 1
    )
    
    if (is.null(n_features)){
        n_features <- nrow(res)
    } else {
        stopifnot(
            length(n_features) == 1,
            round(n_features) == n_features
        )
    }
    cor <- cor(t(assay(res, assay))[, seq_len(n_features)], use = "pairwise.complete.obs")

    range <- max(abs(cor))
    pheatmap(cor, treeheight_row=0, treeheight_col = 0, legend=TRUE,
                    show_colnames=FALSE, show_rownames=TRUE,
                    breaks = seq(-range, range, length.out=100))

}
