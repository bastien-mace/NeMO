#' @include classes.R
#'
#' @title Calculate the Watanabe-Akaike Information Criterion (WAIC)
#'
#' @description This function computes the Watanabe-Akaike Information Criterion
#' (WAIC), also known as the Widely Applicable Information Criterion (Watanabe,
#' 2010). This robust metric is used for model comparison and selection within Bayesian
#' frameworks, offering a way to balance model fit and complexity. The WAIC estimates
#' out-of-sample predictive accuracy, penalising overfitting.
#'
#' @param model A fitted model object created with [Nemodel] with `loglik = TRUE` to record log-likelihoods.
#'
#' @details The WAIC is computed using the following equations, which incorporate both the mean log-likelihood and its variance:
#'
#' - **Log Pointwise Predictive Density (LPPD)**:
#' \deqn{LPPD = \sum_{nijkc=1}^{NIJKC} \ln (\frac{1}{Sim_{tot}} \sum_{sim=1}^{Sim_{tot}} e^{\log \ell_{nijkc,sim}})} \cr
#'   where:
#'   - \eqn{\log \ell_{nijkc,sim}} is the log-likelihood value calculated from simulation \eqn{sim} for species \eqn{n} at site \eqn{i} in sample \eqn{j} in replicate \eqn{k} (except for `'seq_read'` protocol) during campaign. \eqn{c}
#'   - \eqn{Sim_{tot}} is the total number of simulations (MCMC samples) computed.
#' \cr
#' - **Effective Number of Parameters (p_{WAIC})**:
#' \deqn{p_{WAIC} = \sum_{nijkc=1}^{NIJKC} \frac{1}{Sim_{tot}} \sum_{sim=1}^{Sim_{tot}} \left(\log \ell_{nijkc,sim} - \overline{\log \ell_{nijkc}} \right)^2} \cr
#'   where:
#'   - \eqn{\log \ell_{nijkc,sim}} is the log-likelihood value calculated from simulation \eqn{sim} for species \eqn{n} at site \eqn{i} in sample \eqn{j} in replicate \eqn{k} (except for `'seq_read'` protocol) during campaign. \eqn{c}
#'   - \eqn{\overline{\log \ell_{nijkc}}} is the mean log-likelihood value across all simulations for species \eqn{n} at site \eqn{i} in sample \eqn{j} in replicate \eqn{k} (except for `'seq_read'` protocol) during campaign. \eqn{c}
#'   - \eqn{Sim_{tot}} is the total number of simulations (MCMC samples) computed.
#' \cr
#' - **WAIC Score**:
#' \deqn{WAIC = -2 \times LPPD - p_{WAIC}}
#' \cr
#' These computations are performed internally, yielding a scalar WAIC value that summarises the overall performance of the model.
#'
#' @return A structured class object containing three slots:
#' - `waic`: The WAIC score, quantifying the trade-off between model fit and complexity.
#' - `lppd`: The log-predictive density term, reflecting the model's fit.
#' - `p_waic`: The effective number of parameters, representing model complexity.
#'
#' @examples
#' # Load fitted model
#' data(model_PCR_rep)
#'
#' # Calculate the WAIC
#' WAIC(model_PCR_rep)
#'
#' @seealso `Nemodel()`, `covarray()`
#'
#' @export
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats var
#'
WAIC <- function(model)
{

  # Ensure model is provided and has correct type
  if (missing(model)) {
    stop("Error: The 'model' argument is missing.")
  }
  if (!inherits(model, 'Nemodel')) {
    stop("Error: The 'model' argument must be a 'Nemodel' class object.")
  }

  # Ensure protocol is provided and valid
  if (isFALSE(model@loglik)) {
    stop("Error: Log-likelihood was not recorded in the model provided.")
  }


  # Extract log-likelihood values from the model
  loglik <- model@model$BUGSoutput$sims.list$loglik

  iterations <- dim(loglik)[1]

  if(model@protocol == 'PCR_rep'){
    N <- dim(loglik)[2] # Species
    I <- dim(loglik)[3] # Sites
    J <- dim(loglik)[4] # Samples
    K <- dim(loglik)[5] # Replicates
    C <- dim(loglik)[6] # Campaigns

    # Build the log-likelihhod array
    loglik_array <- array(loglik, dim = c(iterations, N, I, J, K, C))

    # Compute the pointwise log-likelihood mean and variance over iterations
    loglik_mean <- apply(exp(loglik_array), c(2, 3, 4, 5, 6), mean)
    loglik_var <- apply(loglik_array, c(2, 3, 4, 5, 6), var)
  }


  if (model@protocol == 'seq_read'){
    I <- dim(loglik)[2] # Sites
    J <- dim(loglik)[3] # Samples
    C <- dim(loglik)[4] # Campaigns

    # Build the log-likelihhod array
    loglik_array <- array(loglik, dim = c(iterations, I, J, C))

    # Compute the pointwise log-likelihood mean and variance over iterations
    loglik_mean <- apply(exp(loglik_array), c(2, 3, 4), mean)
    loglik_var <- apply(loglik_array, c(2, 3, 4), var)
  }


  if (model@protocol == 'PCR_rep_seq_read'){
    I <- dim(loglik)[2] # Sites
    J <- dim(loglik)[3] # Samples
    K <- dim(loglik)[4] # Replicates
    C <- dim(loglik)[5] # Campaigns

    # Build the log-likelihhod array
    loglik_array <- array(loglik, dim = c(iterations, I, J, K, C))

    # Compute the pointwise log-likelihood mean and variance over iterations
    loglik_mean <- apply(exp(loglik_array), c(2, 3, 4, 5), mean)
    loglik_var <- apply(loglik_array, c(2, 3, 4, 5), var)
  }


  # Compute the log pointwise predictive density (lppd)
  lppd <- sum(log(loglik_mean))


  # Compute the effective number of parameters (p_waic)
  p_waic <- sum(loglik_var)


  # Compute WAIC
  waic <- -2 * (lppd - p_waic)


  # Return a class object
  instance <- new('WAIC',
                  waic = waic,
                  lppd = lppd,
                  p_waic = p_waic)

  return(instance)
}
