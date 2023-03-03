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

    expect_true(inherits(p,"pheatmap"))
})
