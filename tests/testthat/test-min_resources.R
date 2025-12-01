test_that("min_ressources works", {
  # Load input array
  data(fish_PCR_rep_seq_read)

  # Load covariate data (Distance to sea)
  data(distance_cov)

  # Build the covariate array applied to sites on the psi level
  covariates <- covarray(protocol = 'PCR_rep_seq_read',
                         array = fish_PCR_rep_seq_read,
                         cov_list = list(Distance = list(cov_data = distance_cov$Distance,
                                                         level = 'psi',
                                                         dimension = 'site')))

  # Run the 'PCR_rep_seq_read' model on the input array with the distance covariate
  model <- Nemodel(protocol = 'PCR_rep_seq_read',
                   array = fish_PCR_rep_seq_read,
                   covariates = covariates)
  unlink('model.txt')

  # Calculate the minimum number of samples, PCR replicates, and sequencing depth
  # required to detect species when present with 95% confidence
  result <- min_resources(model = model,
                          resources = c('J', 'K', 'M'))

  # Test
  expect_equal(all(result@J_min$median == floor(result@J_min$median)), TRUE)
  expect_equal(all(result@J_min$hdi1 == floor(result@J_min$hdi1)), TRUE)
  expect_equal(all(result@J_min$hdi2 == floor(result@J_min$hdi2)), TRUE)
  expect_equal(all(result@K_min$median == floor(result@K_min$median)), TRUE)
  expect_equal(all(result@K_min$hdi1 == floor(result@K_min$hdi1)), TRUE)
  expect_equal(all(result@K_min$hdi2 == floor(result@K_min$hdi2)), TRUE)
  expect_equal(all(result@M_min$median == floor(result@M_min$median)), TRUE)
  expect_equal(all(result@M_min$hdi1 == floor(result@M_min$hdi1)), TRUE)
  expect_equal(all(result@M_min$hdi2 == floor(result@M_min$hdi2)), TRUE)
})
