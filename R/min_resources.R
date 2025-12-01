#' @include classes.R
#'
#' @title Calculate Minimum Resource Requirements for Reliable Species Detection
#'
#' @description This function calculates the minimum resource requirements, *i.e.*
#' samples \eqn{(J_{min})}, PCR replicates \eqn{(K_{min})}, and sequencing depth
#' \eqn{(M_{min})}, needed to confidently detect species when present with a specified
#' confidence level. This estimation relies on established probabilistic models,
#' enabling rigorous resource planning for ecological studies.
#'
#' @param model A fitted [Nemodel] object.
#' @param resources A character vector. Specifies the resource types to compute. Acceptable values include `'J'` (minimum samples), `'K'` (minimum PCR replicates), and `'M'` (minimum sequencing depth). Multiple resource types can be specified simultaneously (*e.g.*, `resources = c('J', 'K', 'M')`).
#' @param conf Numeric. Confidence level, or probability to achieve to detect species when present (default: 0.95).
#'
#' @details The function follows established equations from McArdle (1990) and probabilistic
#' models to compute:
#'
#' - **Minimum number of samples** \eqn{\bold{(J_{min})}}: The function calculates
#' the minimum number of samples required to ensure that the probability of missing
#' species DNA in a sample is below 0.05 with 95% confidence (can be computed for
#' any `protocol`):
#'   \deqn{J_{min} \geq \frac{\ln(0.05)}{\ln(1 - \theta_{nijc})}}{J_min >= ln(0.05) / ln(1 - \theta_nijc)} \cr
#'   where \eqn{\theta_{nijc}} is the probability of DNA collection for species \eqn{n} at site \eqn{i} in sample \eqn{j} during campaign \eqn{c}.
#' \cr
#' - **Minimum number of PCR replicates** \eqn{\bold{(K_{min})}}: The function calculates
#' the minimum number of PCR replicates required to ensure that the probability
#' of missing species DNA in a PCR replicate is below 0.05 with 95% confidence (requires
#' a [Nemodel] object built with `'PCR_rep'` or `'PCR_rep_seq_read'` protocol):
#'   \deqn{K_{min} \geq \frac{\ln(0.05)}{\ln(1 - p_{nijkc})}}{K_min >= ln(0.05) / ln(1 - p_nijkc)} \cr
#'   where \eqn{p_{nijkc}} is the probability of DNA amplification for species \eqn{n} at site \eqn{i} in sample \eqn{j} in replicate \eqn{k} during campaign \eqn{c}.
#' \cr
#' - **Minimum sequencing depth** \eqn{\bold{(M_{min})}}: Following the multinomial
#' theorem, the function calculates the minimum sequencing depth required to ensure
#' that the probability of missing species DNA in a sample or in a PCR replicate
#' is below 0.05 (requires a [Nemodel] object built with `'seq_read'` or `'PCR_rep_seq_read'`
#' protocol):
#'   \deqn{M_{min} \geq \frac{\ln(0.05)}{\ln(1 - \pi_{nijkc})}}{M_min >= ln(0.05) / ln(1 - \pi_nijkc)} \cr
#'   where \eqn{\pi_{nijkc}} is the relative sequence read count for species \eqn{n} at site \eqn{i} in sample \eqn{j} in replicate \eqn{k} (except for `'seq_read'` protocol) during campaign \eqn{c}.
#' \cr
#' The function automatically performs these computations for precise and efficient
#' resource estimation.
#'
#' @return A structured class object with three slots:
#' - `J_min`: Minimum number of samples required to confidently detect species when present across sites.
#' - `K_min`: Minimum number of replicates required to confidently detect species when present across sites and samples.
#' - `M_min`: Minimum sequencing depth required to confidently detect species when present across sites, samples and replicates (except for `'seq_read'` protocol).
#'
#' Each slot contains 3 arrays, corresponding to the median values and the lower and upper bounds of the 95% highest density interval (HDI) of MCMC samples.
#'
#' @examples
#' # Load fitted model
#' data(model_PCR_rep)
#'
#' # Calculate the minimum number of samples and PCR replicates required to
#' # confidently detect species when present with 95% confidence
#' min_resources(model = model_PCR_rep,
#'               resources = c('J', 'K'))
#'
#' @seealso `Nemodel()`, `covarray()`
#'
#' @export
#'
#' @importFrom HDInterval hdi
#' @importFrom magrittr %>%
#' @importFrom methods new
#' @importFrom Rdpack reprompt
#' @importFrom stats median plogis rbinom rnbinom
#'
min_resources <- function(model,
                          resources = c('J'),
                          conf = 0.95)
{

  # Ensure model is provided and has correct type
  if (missing(model)) {
    stop("Error: The 'model' argument is missing.")
  }
  if (!inherits(model, 'Nemodel')) {
    stop("Error: The 'model' argument must be a 'Nemodel' class object.")
  }

  # Define the allowed values for 'resources'
  allowed_resources <- c('J', 'K', 'M')
  if(!is.null(resources) && !all(resources %in% resources)){
    stop("Error: 'resources' can only contain the following values: ",
         paste(allowed_resources, collapse = ', '))
  }

  # Check for numerical arguments
  if (!is.numeric(conf)){
    stop("Error: The confidence level should be numeric.")
  }


  # Set array dimensions
  array <- model@array
  covariates <- model@covariates
  N <- dim(array)[1] # Species
  I <- dim(array)[2] # Sites
  J <- dim(array)[3] # Samples
  if (model@protocol != 'seq_read'){
    K <- dim(array)[4] # Replicates
    C <- dim(array)[5] # Campaigns
  }else{
    C <- dim(array)[4] # Campaigns
  }
  sims <- model@model$BUGSoutput$n.sims # Simulations


  # Initialize minimum resource lists
  J_min = NULL
  K_min = NULL
  M_min = NULL


  # Build the minimum samples list (J_min)
  if ('J' %in% resources){
    # Initialize parameter array
    posterior_theta <- array(0, dim = c(sims, N, I, J, C))
    # Initialize minimum sample lists
    J_min <- list(median = array(NA, dim = c(N, I, J, C)), # Median
                  hdi1 = array(NA, dim = c(N, I, J, C)),   # HDI 1
                  hdi2 = array(NA, dim = c(N, I, J, C)))   # HDI 2

    # Add species covariates estimates to parameter value
    if (!is.null(covariates@theta_cov_sp)){
      Q <- dim(covariates@theta_cov_sp)[2]
      for (q in 1:Q){
        beta_theta_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta_sp[', q, ']'))
        beta_posterior_theta_sp <- model@model$BUGSoutput$sims.matrix[, beta_theta_sp]
        for (n in 1:N){
          for (i in 1:I){
            for (j in 1:J){
              for (c in 1:C){
                posterior_theta[, n, i, j, c] <- posterior_theta[, n, i, j, c] +
                  beta_theta_sp * covariates@theta_cov_sp[n, q]
              }
            }
          }
        }
      }
    }

    # Add other covariates estimates to parameter value
    if (!is.null(covariates@theta_cov)){
      Q <- dim(covariates@theta_cov)[4]
      for (q in 1:Q){
        beta_theta_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta[1,', q, ']'))
        beta_theta_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta[', N, ',', q, ']'))
        beta_posterior_theta <- model@model$BUGSoutput$sims.matrix[, beta_theta_start:beta_theta_stop]
        for (n in 1:N){
          for (i in 1:I){
            for (j in 1:J){
              for (c in 1:C){
                posterior_theta[, n, i, j, c] <- posterior_theta[, n, i, j, c] +
                  beta_posterior_theta[, n] * covariates@theta_cov[i, j, c, q]
              }
            }
          }
        }
      }
    }

    # Add intercept estimate to parameter value
    alpha_theta_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_theta[1]')
    alpha_theta_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_theta[', N, ']'))
    alpha_posterior_theta <- model@model$BUGSoutput$sims.matrix[, alpha_theta_start:alpha_theta_stop]
    for (n in 1:N){
      posterior_theta[, n, , , ] <- posterior_theta[, n, , , ] +
        alpha_posterior_theta[, n]
    }

    # Calculate minimum number of samples
    for (n in 1:N){
      for (i in 1:I){
        for (j in 1:J){
          for (c in 1:C){
            post_theta <- plogis(posterior_theta[, n, i, j, c])
            post_theta <- post_theta[post_theta != 0]
            J_min$median[n, i, j, c] <- ceiling(log(1 - conf) / log(1 - median(post_theta)))
            J_min$hdi1[n, i, j, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_theta)[2]))
            J_min$hdi2[n, i, j, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_theta)[1]))
          }
        }
      }
    }
  }


  # Build the minimum replicates list (K_min)
  if ('K' %in% resources && model@protocol != 'seq_read'){
    # Initialize parameter array
    posterior_p <- array(0, dim = c(sims, N, I, J, K, C))
    # Initialize minimum replicate lists
    K_min <- list(median = array(NA, dim = c(N, I, J, K, C)), # Median
                  hdi1 = array(NA, dim = c(N, I, J, K, C)),   # HDI 1
                  hdi2 = array(NA, dim = c(N, I, J, K, C)))   # HDI 2

    # Add species covariates estimates to parameter value
    if (!is.null(covariates@p_cov_sp)){
      Q <- dim(covariates@p_cov_sp)[2]
      for (q in 1:Q){
        beta_p_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_p_sp[', q, ']'))
        beta_posterior_p_sp <- model@model$BUGSoutput$sims.matrix[, beta_p_sp]
        for (n in 1:N){
          for (i in 1:I){
            for (j in 1:J){
              for (k in 1:K){
                for (c in 1:C){
                  posterior_p[, n, i, j, k, c] <- posterior_p[, n, i, j, k, c] +
                    beta_p_sp * covariates@p_cov_sp[n, q]
                }
              }
            }
          }
        }
      }
    }

    # Add other covariates estimates to parameter value
    if (!is.null(covariates@p_cov)){
      Q <- dim(covariates@p_cov)[5]
      for (q in 1:Q){
        beta_p_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_p[1,', q, ']'))
        beta_p_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_p[', N, ',', q, ']'))
        beta_posterior_p <- model@model$BUGSoutput$sims.matrix[, beta_p_start:beta_p_stop]
        for (n in 1:N){
          for (i in 1:I){
            for (j in 1:J){
              for (k in 1:K){
                for (c in 1:C){
                  posterior_p[, n, i, j, k, c] <- posterior_p[, n, i, j, k, c] +
                    beta_posterior_p[, n] * covariates@p_cov[i, j, k, c, q]
                }
              }
            }
          }
        }
      }
    }

    # Add intercept estimate to parameter value
    alpha_p_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_p[1]')
    alpha_p_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_p[', N, ']'))
    alpha_posterior_p <- model@model$BUGSoutput$sims.matrix[, alpha_p_start:alpha_p_stop]
    for (n in 1:N){
      posterior_p[, n, , , , ] <- posterior_p[, n, , , , ] +
        alpha_posterior_p[, n]
    }

    # Calculate minimum number of replicates
    for (n in 1:N){
      for (i in 1:I){
        for (j in 1:J){
          for (k in 1:K){
            for (c in 1:C){
              post_p <- plogis(posterior_p[, n, i, j, k, c])
              post_p <- post_p[post_p != 0]
              K_min$median[n, i, j, k, c] <- ceiling(log(1 - conf) / log(1 - median(post_p)))
              K_min$hdi1[n, i, j, k, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_p)[2]))
              K_min$hdi2[n, i, j, k, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_p)[1]))
            }
          }
        }
      }
    }
  }


  # Build the minimum sequencing depth list (M_min)
  if ('M' %in% resources && model@protocol != 'PCR_rep'){
    if (model@protocol == 'seq_read'){
      # Build psi posterior distribution: Initialize parameter array
      posterior_psi <- array(0, dim = c(sims, N, I, C))
      # Initialize minimum sequencing depth lists
      M_min <- list(median = array(NA, dim = c(N, I, J, C)), # Median
                    hdi1 = array(NA, dim = c(N, I, J, C)),   # HDI 1
                    hdi2 = array(NA, dim = c(N, I, J, C)))   # HDI 2

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@psi_cov_sp)){
        Q <- dim(covariates@psi_cov_sp)[2]
        for (q in 1:Q){
          beta_psi_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_psi_sp[', q, ']'))
          beta_posterior_psi_sp <- model@model$BUGSoutput$sims.matrix[, beta_psi_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (c in 1:C){
                posterior_psi[, n, i, c] <- posterior_psi[, n, i, c] +
                  beta_psi_sp * covariates@psi_cov_sp[n, q]
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@psi_cov)){
        Q <- dim(covariates@psi_cov)[3]
        for (q in 1:Q){
          beta_psi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_psi[1,', q, ']'))
          beta_psi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_psi[', N, ',', q, ']'))
          beta_posterior_psi <- model@model$BUGSoutput$sims.matrix[, beta_psi_start:beta_psi_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (c in 1:C){
                posterior_psi[, n, i, c] <- posterior_psi[, n, i, c] +
                  beta_posterior_psi[, n] * covariates@psi_cov[i, c, q]
              }
            }
          }
        }
      }

      # Add intercept estimate to parameter value
      alpha_psi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_psi[1]')
      alpha_psi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_psi[', N, ']'))
      alpha_posterior_psi <- model@model$BUGSoutput$sims.matrix[, alpha_psi_start:alpha_psi_stop]
      for (n in 1:N){
        posterior_psi[, n, ,] <- posterior_psi[, n, ,] +
          alpha_posterior_psi[, n]
      }


      # Build theta posterior distribution: Initialize parameter array
      posterior_theta <- array(0, dim = c(sims, N, I, J, C))

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@theta_cov_sp)){
        Q <- dim(covariates@theta_cov_sp)[2]
        for (q in 1:Q){
          beta_theta_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta_sp[', q, ']'))
          beta_posterior_theta_sp <- model@model$BUGSoutput$sims.matrix[, beta_theta_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (c in 1:C){
                  posterior_theta[, n, i, j, c] <- posterior_theta[, n, i, j, c] +
                    beta_theta_sp * covariates@theta_cov_sp[n, q]
                }
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@theta_cov)){
        Q <- dim(covariates@theta_cov)[4]
        for (q in 1:Q){
          beta_theta_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta[1,', q, ']'))
          beta_theta_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta[', N, ',', q, ']'))
          beta_posterior_theta <- model@model$BUGSoutput$sims.matrix[, beta_theta_start:beta_theta_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (c in 1:C){
                  posterior_theta[, n, i, j, c] <- posterior_theta[, n, i, j, c] +
                    beta_posterior_theta[, n] * covariates@theta_cov[i, j, c, q]
                }
              }
            }
          }
        }
      }

      # Add intercept estimate to parameter value
      alpha_theta_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_theta[1]')
      alpha_theta_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_theta[', N, ']'))
      alpha_posterior_theta <- model@model$BUGSoutput$sims.matrix[, alpha_theta_start:alpha_theta_stop]
      for (n in 1:N){
        posterior_theta[, n, , ,] <- posterior_theta[, n, , ,] +
          alpha_posterior_theta[, n]
      }

      # Build phi posterior distribution: Initialize parameter array
      posterior_phi <- array(0, dim = c(sims, N, I, J, C))

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@phi_cov_sp)){
        Q <- dim(covariates@phi_cov_sp)[2]
        for (q in 1:Q){
          beta_phi_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_phi_sp[', q, ']'))
          beta_posterior_phi_sp <- model@model$BUGSoutput$sims.matrix[, beta_phi_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (c in 1:C){
                  posterior_phi[, n, i, j, c] <- posterior_phi[, n, i, j, c] +
                    beta_phi_sp * covariates@phi_cov_sp[n, q]
                }
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@phi_cov)){
        Q <- dim(covariates@phi_cov)[4]
        for (q in 1:Q){
          beta_phi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_phi[1,', q, ']'))
          beta_phi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_phi[', N, ',', q, ']'))
          beta_posterior_phi <- model@model$BUGSoutput$sims.matrix[, beta_phi_start:beta_phi_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (c in 1:C){
                  posterior_phi[, n, i, j, c] <- posterior_phi[, n, i, j, c] +
                    beta_posterior_phi[, n] * covariates@phi_cov[i, j, c, q]
                }
              }
            }
          }
        }
      }

      # Add intercept estimate to parameter value
      alpha_phi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_phi[1]')
      alpha_phi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_phi[', N, ']'))
      alpha_posterior_phi <- model@model$BUGSoutput$sims.matrix[, alpha_phi_start:alpha_phi_stop]
      for (n in 1:N){
        posterior_phi[, n, , ,] <- posterior_phi[, n, , ,] +
          alpha_posterior_phi[, n]
      }


      # Initialize Z, A, and S arrays
      Z <- array(0, dim = c(sims, N, I, C))
      A <- array(0, dim = c(sims, N, I, J, C))
      S <- array(0, dim = c(sims, N, I, J, C))

      # Initialize pi array
      posterior_pi <- array(0, dim = c(sims, N, I, J, C))
      delta <- model@model$BUGSoutput$sims.matrix[, 'delta']
      for (n in 1:N){
        for (c in 1:C){
          for (i in 1:I){
            # Rebuild Z array
            Z[, n, i, c] <- rbinom(n = sims, size = 1, prob = plogis(posterior_psi[, n, i, c]))
            for(j in 1:J){
              # Rebuild A array
              A[, n, i, j, c] <- rbinom(n = sims, size = 1, prob = Z[, n, i, c] * plogis(posterior_theta[, n, i, j, c]))
              # Rebuild S array
              S[, n, i, j, c] <- rnbinom(n = sims, size = delta, prob = delta / (delta + exp(posterior_phi[, n, i, j, c])))
              # Rebuild pi distribution
              posterior_pi[, n, i, j, c] <- A[, n, i, j, c] * S[, n, i, j, c] / ifelse(apply(A[, , i, j, c] * S[, , i, j, c], c(1), sum, na.rm = T) > 0,
                                                                                       apply(A[, , i, j, c] * S[, , i, j, c], c(1), sum, na.rm = T),
                                                                                       1)
              post_pi <- posterior_pi[, n, i, j, c]
              # Remove null values
              post_pi <- post_pi[post_pi != 0]

              # Calculate minimum sequencing depth
              M_min$median[n, i, j, c] <- ceiling(log(1 - conf) / log(1 - median(post_pi)))
              M_min$hdi1[n, i, j, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_pi)[2]))
              M_min$hdi2[n, i, j, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_pi)[1]))
            }
          }
        }
      }
    }


    if (model@protocol == 'PCR_rep_seq_read'){
      # Build psi posterior distribution: Initialize parameter array
      posterior_psi <- array(0, dim = c(sims, N, I, C))
      # Initialize minimum sequencing depth lists
      M_min <- list(median = array(NA, dim = c(N, I, J, K, C)), # Median
                    hdi1 = array(NA, dim = c(N, I, J, K, C)),   # HDI 1
                    hdi2 = array(NA, dim = c(N, I, J, K, C)))   # HDI 2

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@psi_cov_sp)){
        Q <- dim(covariates@psi_cov_sp)[2]
        for (q in 1:Q){
          beta_psi_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_psi_sp[', q, ']'))
          beta_posterior_psi_sp <- model@model$BUGSoutput$sims.matrix[, beta_psi_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (c in 1:C){
                posterior_psi[, n, i, c] <- posterior_psi[, n, i, c] +
                  beta_psi_sp * covariates@psi_cov_sp[n, q]
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@psi_cov)){
        Q <- dim(covariates@psi_cov)[3]
        for (q in 1:Q){
          beta_psi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_psi[1,', q, ']'))
          beta_psi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_psi[', N, ',', q, ']'))
          beta_posterior_psi <- model@model$BUGSoutput$sims.matrix[, beta_psi_start:beta_psi_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (c in 1:C){
                posterior_psi[, n, i, c] <- posterior_psi[, n, i, c] +
                  beta_posterior_psi[, n] * covariates@psi_cov[i, c, q]
              }
            }
          }
        }
      }

      # Add intercept estimate to parameter value
      alpha_psi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_psi[1]')
      alpha_psi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_psi[', N, ']'))
      alpha_posterior_psi <- model@model$BUGSoutput$sims.matrix[, alpha_psi_start:alpha_psi_stop]

      for (n in 1:N){
        posterior_psi[, n, ,] <- posterior_psi[, n, ,] +
          alpha_posterior_psi[, n]
      }


      # Build theta posterior distribution: Initialize parameter array
      posterior_theta <- array(0, dim = c(sims, N, I, J, C))

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@theta_cov_sp)){
        Q <- dim(covariates@theta_cov_sp)[2]
        for (q in 1:Q){
          beta_theta_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta_sp[', q, ']'))
          beta_posterior_theta_sp <- model@model$BUGSoutput$sims.matrix[, beta_theta_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (c in 1:C){
                  posterior_theta[, n, i, j, c] <- posterior_theta[, n, i, j, c] +
                    beta_theta_sp * covariates@theta_cov_sp[n, q]
                }
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@theta_cov)){
        Q <- dim(covariates@theta_cov)[4]
        for (q in 1:Q){
          beta_theta_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta[1,', q, ']'))
          beta_theta_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_theta[', N, ',', q, ']'))
          beta_posterior_theta <- model@model$BUGSoutput$sims.matrix[, beta_theta_start:beta_theta_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (c in 1:C){
                  posterior_theta[, n, i, j, c] <- posterior_theta[, n, i, j, c] +
                    beta_posterior_theta[, n] * covariates@theta_cov[i, j, c, q]
                }
              }
            }
          }
        }
      }

      # Add intercept estimate to parameter value
      alpha_theta_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_theta[1]')
      alpha_theta_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_theta[', N, ']'))
      alpha_posterior_theta <- model@model$BUGSoutput$sims.matrix[, alpha_theta_start:alpha_theta_stop]
      for (n in 1:N){
        posterior_theta[, n, , ,] <- posterior_theta[, n, , ,] +
          alpha_posterior_theta[, n]
      }


      # Build p posterior distribution: Initialize parameter array
      posterior_p <- array(0, dim = c(sims, N, I, J, K, C))

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@p_cov_sp)){
        Q <- dim(covariates@p_cov_sp)[2]
        for (q in 1:Q){
          beta_p_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_p_sp[', q, ']'))
          beta_posterior_p_sp <- model@model$BUGSoutput$sims.matrix[, beta_p_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (k in 1:K){
                  for (c in 1:C){
                    posterior_p[, n, i, j, k, c] <- posterior_p[, n, i, j, k, c] +
                      beta_p_sp * covariates@p_cov_sp[n, q]
                  }
                }
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@p_cov)){
        Q <- dim(covariates@p_cov)[5]
        for (q in 1:Q){
          beta_p_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_p[1,', q, ']'))
          beta_p_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_p[', N, ',', q, ']'))
          beta_posterior_p <- model@model$BUGSoutput$sims.matrix[, beta_p_start:beta_p_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (k in 1:K){
                  for (c in 1:C){
                    posterior_p[, n, i, j, k, c] <- posterior_p[, n, i, j, k, c] +
                      beta_posterior_p[, n] * covariates@p_cov[i, j, k, c, q]
                  }
                }
              }
            }
          }
        }
      }

      # Add intercept estimate to parameter value
      alpha_p_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_p[1]')
      alpha_p_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_p[', N, ']'))
      alpha_posterior_p <- model@model$BUGSoutput$sims.matrix[, alpha_p_start:alpha_p_stop]
      for (n in 1:N){
        posterior_p[, n, , , ,] <- posterior_p[, n, , , ,] +
          alpha_posterior_p[, n]
      }

      # Build phi posterior distribution: Initialize parameter array
      posterior_phi <- array(0, dim = c(sims, N, I, J, K, C))

      # Add species covariates estimates to parameter value
      if (!is.null(covariates@phi_cov_sp)){
        Q <- dim(covariates@phi_cov_sp)[2]
        for (q in 1:Q){
          beta_phi_sp <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_phi_sp[', q, ']'))
          beta_posterior_phi_sp <- model@model$BUGSoutput$sims.matrix[, beta_phi_sp]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (k in 1:K){
                  for (c in 1:C){
                    posterior_phi[, n, i, j, k, c] <- posterior_phi[, n, i, j, k, c] +
                      beta_phi_sp * covariates@phi_cov_sp[n, q]
                  }
                }
              }
            }
          }
        }
      }

      # Add other covariates estimates to parameter value
      if (!is.null(covariates@phi_cov)){
        Q <- dim(covariates@phi_cov)[5]
        for (q in 1:Q){
          beta_phi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_phi[1,', q, ']'))
          beta_phi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('beta_phi[', N, ',', q, ']'))
          beta_posterior_phi <- model@model$BUGSoutput$sims.matrix[, beta_phi_start:beta_phi_stop]
          for (n in 1:N){
            for (i in 1:I){
              for (j in 1:J){
                for (k in 1:K){
                  for (c in 1:C){
                    posterior_phi[, n, i, j, k, c] <- posterior_phi[, n, i, j, k, c] +
                      beta_posterior_phi[, n] * covariates@phi_cov[i, j, k, c, q]
                  }
                }
              }
            }
          }
        }
      }

      # Add parameter estimate to parameter value
      alpha_phi_start <- which(colnames(model@model$BUGSoutput$sims.matrix) == 'alpha_phi[1]')
      alpha_phi_stop <- which(colnames(model@model$BUGSoutput$sims.matrix) == paste0('alpha_phi[', N, ']'))
      alpha_posterior_phi <- model@model$BUGSoutput$sims.matrix[, alpha_phi_start:alpha_phi_stop]
      for (n in 1:N){
        posterior_phi[, n, , , ,] <- posterior_phi[, n, , , ,] +
          alpha_posterior_phi[, n]
      }


      # Initialize Z, A, W, and S arrays
      Z <- array(0, dim = c(sims, N, I, C))
      A <- array(0, dim = c(sims, N, I, J, C))
      W <- array(0, dim = c(sims, N, I, J, K, C))
      S <- array(0, dim = c(sims, N, I, J, K, C))

      # Initialize pi array
      posterior_pi <- array(0, dim = c(sims, N, I, J, K, C))
      delta <- model@model$BUGSoutput$sims.matrix[, 'delta']
      for (n in 1:N){
        for (c in 1:C){
          for (i in 1:I){
            # Rebuild Z array
            Z[, n, i, c] <- rbinom(n = sims, size = 1, prob = plogis(posterior_psi[, n, i, c]))
            for(j in 1:J){
              # Rebuild A array
              A[, n, i, j, c] <- rbinom(n = sims, size = 1, prob = Z[, n, i, c] * plogis(posterior_theta[, n, i, j, c]))
              for (k in 1:K){
                # Rebuild W array
                W[, n, i, j, k, c] <- rbinom(n = sims, size = 1,prob = A[, n, i, j, c] * plogis(posterior_p[, n, i, j, k, c]))
                # Rebuild S array
                S[, n, i, j, k, c] <- rnbinom(n = sims, size = delta, prob = delta / (delta + exp(posterior_phi[, n, i, j, k, c])))
                # Rebuild pi distribution
                posterior_pi[, n, i, j, k, c] <- W[, n, i, j, k, c] * S[, n, i, j, k, c] / ifelse(apply(W[, , i, j, k, c] * S[, , i, j, k, c], c(1), sum, na.rm = T) > 0,
                                                                                                  apply(W[, , i, j, k, c] * S[, , i, j, k, c], c(1), sum, na.rm = T),
                                                                                                  1)
                post_pi <- posterior_pi[, n, i, j, k, c]
                # Remove null values
                post_pi <- post_pi[post_pi != 0]

                # Calculate minimum sequencing depth
                M_min$median[n, i, j, k, c] <- ceiling(log(1 - conf) / log(1 - median(post_pi)))
                M_min$hdi1[n, i, j, k, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_pi)[2]))
                M_min$hdi2[n, i, j, k, c] <- ceiling(log(1 - conf) / log(1 - hdi(post_pi)[1]))
              }
            }
          }
        }
      }
    }
  }

  # Return a class object
  instance <- new('min_resources',
                  J_min = J_min,
                  K_min = K_min,
                  M_min = M_min)

  return(instance)
}
