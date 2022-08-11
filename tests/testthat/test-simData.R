test_that("simData returns correct output",{
  samples <- 100
  latent_factors <- 3
  replicates <- 2

  data <- simData(samples=samples,latent_factors =latent_factors,
                  replicates = replicates)

  expect_equal(dim(data)[1], samples)
  expect_equal(dim(data)[2], latent_factors*10)
  expect_equal(dim(data)[3], replicates)
  expect_true(inherits(data, "array"))
})
