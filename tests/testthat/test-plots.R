test_that("plotSelectedFeatures returns plot",{
    data <-  array(rnorm(50*30), dim=c(30,50),
                dimnames=list(paste("feature", seq_len(30)), NULL))
    reps <- rep(c(1,2), each=25)

    selected <- paste("feature", seq_len(5))
    res <- data.frame(selected = selected)

    se_res <- SummarizedExperiment::SummarizedExperiment(
        assays=list(selected=data[selected,]),
        rowData=res)

    p <- plotSelectedFeatures(se_res)
    
    expect_error(plotSelectedFeatures(data))
    expect_error(plotSelectedFeatures(se_res, n_features=c("feature 1", "feature 2")))
    expect_error(plotSelectedFeatures(se_res, assay="sel_assay"))
    
    expect_true(all(p$tree_row$labels == selected))
    expect_true(is(p,"pheatmap"))
})

test_that("plotVarianceExplained returns plot",{
    data <-  array(rnorm(50*30), dim=c(30,50),
                   dimnames=list(paste("feature", seq_len(30)), NULL))
    reps <- rep(c(1,2), each=25)
    
    selected <- paste("feature", seq_len(5))
    res <- data.frame(selected = selected)
    se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(selected=data[selected,]),
      rowData=res)
    
    se_res <- FeatSeek(data, reps, max_features=5)
   
    expect_true(is.null(plotVarianceExplained(se_res)))
    expect_error(plotVarianceExplained(data))
    expect_error(plotVarianceExplained(se), "'res' does not contain attribute explained_variance.  Please run FeatSeek first!")
    
    
})
  