test_that("plotSelectedFeatures returns plot",{
    data <-  array(rnorm(50*30*2), dim=c(50,30,2), dimnames=list(NULL, paste("feature", seq_len(30)), NULL))
    selected <- paste("feature", seq_len(5))
    res <- data.frame(selected = selected)
    p <- plotSelectedFeatures(data, res)

    expect_true(inherits(p,"pheatmap"))
})
