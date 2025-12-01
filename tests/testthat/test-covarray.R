test_that("covtable works", {
  # Load input array
  data(fish_PCR_rep)

  # Load covariate data (Distance to sea)
  data(distance_cov)

  # Build the covariate array applied to sites on the psi level
  result <- covarray(protocol = 'PCR_rep',
                     array = fish_PCR_rep,
                     cov_list = list(Distance = list(cov_data = distance_cov$Distance,
                                                     level = 'psi',
                                                     dimension = 'site')))

  # Test
  expect_equal(all(is.numeric(result@psi_cov)), TRUE)
})
