test_that("WAIC works", {
  # Load input array
  data(fish_PCR_rep)

  # Load covariate data (Distance to sea)
  data(distance_cov)

  # Build the covariate array applied to sites on the psi level
  covariates <- covarray(protocol = 'PCR_rep',
                         array = fish_PCR_rep,
                         cov_list = list(Distance = list(cov_data = distance_cov$Distance,
                                                         level = 'psi',
                                                         dimension = 'site')))

  # Run the 'PCR_rep' model on the input array with the distance covariate,
  # and record the log-likelihoods for model comparison
  model <- Nemodel(protocol = 'PCR_rep',
                   array = fish_PCR_rep,
                   covariates = covariates,
                   loglik = TRUE)
  unlink('model.txt')

  # Calculate the WAIC
  result <- WAIC(model)

  # Test
  expect_equal(is.numeric(result@waic), TRUE)
  expect_equal(is.numeric(result@lppd), TRUE)
  expect_equal(is.numeric(result@p_waic), TRUE)
})
