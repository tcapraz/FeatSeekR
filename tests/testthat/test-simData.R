test_that("simData returns correct output",{
    samples <- 100
    latent_factors <- 3
    replicates <- 2

    sim <- simData(samples=samples,latent_factors =latent_factors,
                    replicates = replicates)
    data <- sim[[1]]
    reps <- sim[[2]]
    expect_equal(dim(data)[1], samples*replicates)
    expect_equal(dim(data)[2], latent_factors*10)
    expect_equal(length(reps), replicates*samples)
    expect_true(inherits(data, "array"))
})
