test_that("simData returns correct output",{
    conditions <- 100
    latent_factors <- 3
    replicates <- 2

    sim <- simData(conditions=conditions,n_latent_factors =latent_factors,
                    replicates = replicates)
    data <- SummarizedExperiment::assays(sim)$data
    conds <- SummarizedExperiment::colData(sim)$conditions
    
    expect_error(simData(conditions=conditions, n_latent_factors=latent_factors,
                         replicates="three"))
    expect_error(simData(conditions=conditions, n_latent_factors=latent_factors,
                         replicates=c(1,2,3)))
    expect_error(simData(conditions=conditions, n_latent_factors=c(1,2,3),
                         replicates=conds))
    
    expect_equal(dim(data)[2], conditions*replicates)
    expect_equal(dim(data)[1], latent_factors*10)
    expect_equal(length(conds), replicates*conditions)
    expect_true(all(conds == rep(seq_len(conditions), replicates)))
    expect_true(is(sim, "SummarizedExperiment"))
})
