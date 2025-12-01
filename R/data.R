#' Sequence Read Count Data with Individually Indexed PCR Replicates
#'
#' A 5D array containing sequence read count data for fish species across sites, samples, PCR replicates, and campaigns.
#'
#' @format A 5-dimensional array with:
#' \describe{
#'   \item{Dimension 1}{Species (10)}
#'   \item{Dimension 2}{Sites (10)}
#'   \item{Dimension 3}{Samples (2)}
#'   \item{Dimension 4}{PCR Replicates (5)}
#'   \item{Dimension 5}{Campaigns (1)}
#' }
#' This dataset can directly be used in [Nemodel] with the `'PCR_rep_seq_read'` protocol.
#' @usage data(fish_PCR_rep_seq_read)
#' @seealso \code{\link{Nemodel}}
'fish_PCR_rep_seq_read'

#' Presence/Absence Data with Individually Indexed PCR Replicates
#'
#' A 5D array containing presence/absence data for fish species across sites, samples, PCR replicates, and campaigns.
#'
#' @format A 5-dimensional array with:
#' \describe{
#'   \item{Dimension 1}{Species (10)}
#'   \item{Dimension 2}{Sites (10)}
#'   \item{Dimension 3}{Samples (2)}
#'   \item{Dimension 4}{PCR Replicates (5)}
#'   \item{Dimension 5}{Campaigns (1)}
#' }
#' This dataset can directly be used in [Nemodel] with the `'PCR_rep'` protocol.
#' @usage data(fish_PCR_rep)
#' @seealso \code{\link{Nemodel}}
'fish_PCR_rep'

#' Sequence Read Count Data with Pooled PCR Replicates
#'
#' A 4D array containing sequence read count data for fish species across sites, samples, and campaigns.
#'
#' @format A 4-dimensional array with:
#' \describe{
#'   \item{Dimension 1}{Species (10)}
#'   \item{Dimension 2}{Sites (10)}
#'   \item{Dimension 3}{Samples (2)}
#'   \item{Dimension 4}{Campaigns (1)}
#' }
#' This dataset can directly be used in [Nemodel] with the `'seq_read'` protocol.
#' @usage data(fish_seq_read)
#' @seealso \code{\link{Nemodel}}
'fish_seq_read'

#' Distance to Sea Covariate Data
#'
#' A dataframe containing the standardised distance to the sea for each of the 10 sites in the fish datasets.
#'
#' @format A dataframe with 10 rows and 1 column:
#' \describe{
#'   \item{Distance}{Standardised distance to the sea for each site}
#' }
#' The `$Distance` column can directly be used in [covarray].
#' @usage data(distance_cov)
#' @seealso \code{\link{covarray}}
'distance_cov'

#' Example model output using the `PCR_rep` protocol
#'
#' This dataset contains the output from running [Nemodel] on [fish_PCR_rep] integrating [distance_cov] as covariate through [covarray], using the `PCR_rep` protocol and the following parameters: `nb_iterations = 300`, `nb_burnin = 150`, `nb_thinning = 1`, `nb_chains = 2`.
#'
#' @format A structured class object ([Nemodel] object) with six slots:
#' \describe{
#'   \item{model}{A fitted model in `rjags` format}
#'   \item{protocol}{The modelling protocol used (`PCR_rep`)}
#'   \item{array}{The model input data array (`fish_PCR_rep`)}
#'   \item{covariates}{The covariate list implemented in the model (built from `distance_cov`)}
#'   \item{names}{A list containing species, site, sample, replicate, campaign and estimates' names}
#'   \item{loglik}{\code{TRUE}: log-likelihood was stored for model comparison}
#' }
#'
#' @usage data(model_PCR_rep)
#' @seealso \code{\link{Nemodel}}, \code{\link{covarray}}, \code{\link{fish_PCR_rep}}, \code{\link{distance_cov}}
"model_PCR_rep"
