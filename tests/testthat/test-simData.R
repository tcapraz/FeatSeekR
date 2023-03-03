test_that("simData returns correct output",{
    samples <- 100
    latent_factors <- 3
    replicates <- 2

    sim <- simData(samples=samples,n_latent_factors =latent_factors,
                    replicates = replicates)
    data <- SummarizedExperiment::assays(sim)$data
    reps <- SummarizedExperiment::colData(sim)$replicates
    expect_equal(dim(data)[2], samples*replicates)
    expect_equal(dim(data)[1], latent_factors*10)
    expect_equal(length(reps), replicates*samples)
    expect_true(inherits(sim, "SummarizedExperiment"))
})
