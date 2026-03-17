test_that("Nemodel works", {
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

  # Run the 'PCR_rep' model on the input array with the distance covariate
  result <- capture.output(Nemodel(protocol = 'PCR_rep',
                                   array = fish_PCR_rep,
                                   covariates = covariates))
  unlink('model.txt')

  # Test
  var <- strsplit(result, "\n")
  expect_equal("An object of class \"Nemodel\"" %in% var, TRUE)
})

