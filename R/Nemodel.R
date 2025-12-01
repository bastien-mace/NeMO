#' @include classes.R
#'
#' @title Bayesian Multispecies Occupancy Model for eDNA Data.
#'
#' @description This function is the core of the **NeMO** package, designed for
#' Bayesian modelling of multispecies occupancy in eDNA metabarcoding studies. It
#' provides flexibility to fit different modelling protocols, accommodating various
#' study designs and data structures. The function uses the **R2jags** package for
#' MCMC sampling.
#'
#' @param protocol Character string. Specifies the modeling protocol. \cr
#'   Options are:
#'   - `'PCR_rep'`: Requires a 5D (N × I × J × K × C) presence/absence input array
#'   - `'seq_read'`: Requires a 4D (N × I × J × C) sequence read count input array
#'   - `'PCR_rep_seq_read'`: Requires a 5D (N × I × J × K × C) sequence read count input array \cr
#'     where:
#'      - N: Number of species.
#'      - I: Number of sites.
#'      - J: Number of samples.
#'      - K: Number of PCR replicates.
#'      - C: Number of campaigns.
#' @param array Input data array. Its required dimensionality depends on the selected `protocol` argument.
#' @param covariates Optional. A preprocessed covariate list generated using [covarray].
#' @param name Character string specifying the filename and path for recording the BUGS language model file.
#' @param nb_iterations Integer. Total number of MCMC iterations to run (default: 2000).
#' @param nb_burnin Integer. Number of initial burn-in iterations (default: `floor(nb_iterations / 2)`).
#' @param nb_thinning Integer. Thinning interval to reduce autocorrelation in MCMC samples (default: `max(1, floor((nb_iterations - nb_burnin) / 1000))`).
#' @param nb_chains Integer. Number of independent MCMC chains to compute (default: 3).
#' @param parallel Logical. If `TRUE`, enables parallel computation to speed up sampling for large datasets (default: `FALSE`).
#' @param loglik Logical. If `TRUE`, computes and saves log-likelihoods for model comparison using [WAIC] (default: `FALSE`).
#' @param latent Optional character vector. Specifies latent arrays (*e.g.*, `'Z'`, `'A'`, `'W'`, `'S'`, `'Y'`) to save.
#' @param posterior Optional character vector. Specifies posterior distributions of key parameters (*e.g.*, `'psi'`, `'theta'`, `'p'`, `'pi'`, `'phi'`) to save.
#' @param tau Numeric. Precision parameter for normal (hyper)priors (default: 1).
#' @param rho Numeric. Shape parameter of the Gamma prior (only for `'seq_read'` or `'PCR_rep_seq_read'` protocols, default: 1).
#' @param lambda Numeric. Rate parameter of the Gamma prior (only for `'seq_read'` or `'PCR_rep_seq_read'` protocols, default: 1).
#' @param size Numeric. Initialization size parameter for the \eqn{S} array (only for `'seq_read'` or `'PCR_rep_seq_read'` protocols, default: 1).
#' @param prob Numeric. Initialization probability parameter for the \eqn{S} array (only for `'seq_read'` or `'PCR_rep_seq_read'` protocols, default: 0.001).
#' @param ... Additional arguments (related to [jags] or [jags.parallel] functions from **R2jags**).
#'
#' @return A structured class object ([Nemodel] object) with five slots:
#'   - `model`: The fitted model output in `rjags` format, with a BUGS language text file stored at the specified location (`name` argument).
#'   - `protocol`: The modelling protocol used (*i.e.*, `'PCR_rep'`, `'seq_read'`, or `'PCR_rep_seq_read'`).
#'   - `array`: The model input data array.
#'   - `covariates`: The covariate list implemented in the model.
#'   - `names`: A list containing species, site, sample, replicate, and campaign names. This slot also contains:
#'       - `estimates_sp`: Names associated with the species-level random effects.
#'       - `estimates`: Names associated with the intercepts and other random effects.
#'   - `loglik`: Logical. If \code{TRUE}, log-likelihood values are stored for use in model comparison.
#'
#' @details This function provides flexibility for different study designs by allowing users to model occupancy using either presence/absence arrays (`PCR_rep`) or sequence read counts (`seq_read`, `PCR_rep_seq_read`).
#' Covariates can be incorporated to account for external predictors.
#' The MCMC sampling process can be customized using multiple control parameters (`nb_iterations`, `nb_burnin`, `nb_thinning`, etc.), and parallel computation can be enabled with `parallel = TRUE`.
#' The function also allows users to store log-likelihoods, latent arrays, and posterior distributions for downstream analysis.
#'
#' @examples
#' # Load input array
#' data(fish_PCR_rep)
#'
#' # Load covariate data (Distance to sea)
#' data(distance_cov)
#'
#' # Build the covariate array applied to sites on the psi level
#' covariates <- covarray(protocol = 'PCR_rep',
#'                        array = fish_PCR_rep,
#'                        cov_list = list(Distance = list(cov_data = distance_cov$Distance,
#'                                                        level = 'psi',
#'                                                        dimension = 'site')))
#'
#' # Run the 'PCR_rep' model on the input array with the distance covariate
#' Nemodel(protocol = 'PCR_rep',
#'         array = fish_PCR_rep,
#'         covariates = covariates)
#' unlink('model.txt')
#'
#' @seealso \code{\link{WAIC}}, \code{\link{covarray}}
#'
#' @export
#'
#' @importFrom methods new
#' @importFrom R2jags jags jags.parallel
#' @importFrom Rdpack reprompt
#' @importFrom stats rnbinom
#' @importFrom utils tail
#'
Nemodel <- function(protocol = c('PCR_rep', 'seq_read', 'PCR_rep_seq_read'),
                    array,
                    covariates = NULL,
                    name = 'model',
                    nb_iterations = 1000,
                    nb_burnin = floor(nb_iterations / 2),
                    nb_thinning = max(1, floor((nb_iterations - nb_burnin) / 1000)),
                    nb_chains = 2,
                    parallel = FALSE,
                    loglik = FALSE,
                    latent = NULL,
                    posterior = NULL,
                    tau = 1,
                    rho = 1,
                    lambda = 1,
                    size = 1,
                    prob = 1e-3, ...)
{

 # Ensure protocol is provided and valid
  if (missing(protocol)) {
    stop("Error: The 'protocol' argument is missing.")
  }
  protocol <- match.arg(protocol)

  # Ensure array is provided and has correct dimensions
  if (missing(array)) {
    stop("Error: The 'array' argument is missing.")
  }

  if (!is.array(array)) {
    stop("Error: 'array' must be an array with the required dimensions.")
  }

  if (protocol == 'PCR_rep' && length(dim(array)) != 5) {
    stop("Error: 'array' must have 5 dimensions (N, I, J, K, C) for 'PCR_rep' protocol.")
  }
  if (protocol == 'seq_read' && length(dim(array)) != 4) {
    stop("Error: 'array' must have 4 dimensions (N, I, J, C) for 'seq_read' protocol.")
  }
  if (protocol == 'PCR_rep_seq_read' && length(dim(array)) != 5) {
    stop("Error: 'array' must have 5 dimensions (N, I, J, K, C) for 'PCR_rep_seq_read' protocol.")
  }

  # Check for correct covariates format
  if (!is.null(covariates) && !inherits(covariates, 'covarray')) {
    stop("Error: The 'covariates' argument must be a 'covarray' class object.")
  }

  # Check latent variables
  allowed_latents <- c('Z', 'A', 'W', 'S', 'Y')
  if (!is.null(latent) && !all(latent %in% allowed_latents)) {
    stop("Error: 'latent' can only contain the following values: ",
         paste(allowed_latents, collapse = ', '))
  }

  # Check posterior variables
  allowed_posteriors <- c('psi', 'theta', 'p', 'phi', 'pi')
  if (!is.null(posterior) && !all(posterior %in% allowed_posteriors)) {
    stop("Error: 'posterior' can only contain the following values: ",
         paste(allowed_posteriors, collapse = ', '))
  }

  # Check for dimension names
  if (is.null(dimnames(array)[[1]])) {
    stop("Error: Please provide names for dimension 1 (species).")
  }
  if (is.null(dimnames(array)[[2]])) {
    stop("Error: Please provide names for dimension 2 (sites).")
  }
  if (is.null(dimnames(array)[[3]])) {
    stop("Error: Please provide names for dimension 3 (samples).")
  }
  if (protocol != 'seq_read') {
    if (is.null(dimnames(array)[[4]])) {
      stop("Error: Please provide names for dimension 4 (replicates).")
    }
    if (is.null(dimnames(array)[[5]])) {
      stop("Error: Please provide names for dimension 5 (campaigns).")
    }
  } else {
    if (is.null(dimnames(array)[[4]])) {
      stop("Error: Please provide names for dimension 4 (campaigns).")
    }
  }

  # Check for numerical arguments
  if (tau <= 0){
    stop("Error: 'tau' should be positive.")
  }
  if (rho <= 0){
    stop("Error: 'rho' should be positive.")
  }
  if (lambda <= 0){
    stop("Error: 'lambda' should be positive.")
  }
  if (size <= 0){
    stop("Error: 'size' should be positive.")
  }
  if (prob <= 0 | prob > 1){
    stop("Error: 'prob' should be in ]0,1].")
  }


  # Set input array dimensions
  Y <- array
  N <- dim(Y)[1]   # Species
  I <- dim(Y)[2]   # Sites
  J <- dim(Y)[3]   # Samples
  if (protocol != 'seq_read'){
    K <- dim(Y)[4] # Replicates
    C <- dim(Y)[5] # Campaigns
  }else{
    C <- dim(Y)[4] # Campaigns
  }


  # Calculate sequencing depth and handle NA values
  if(protocol == 'seq_read'){
    M <- apply(Y, c(2, 3, 4), sum, na.rm = T)    # Sequence read count per sample
                                                 # (proxy for sequencing depth)
    # Set Y to NA when sequencing depth is 0 (sequencing failed for all species)
    for (i in 1:I){
      for (j in 1:J){
        for (c in 1:C){
          if (M[i, j, c] == 0){
            Y[, i, j, c] <- NA
          }
        }
      }
    }
  }

  if (protocol == 'PCR_rep_seq_read'){
    M <- apply(Y, c(2, 3, 4, 5), sum, na.rm = T) # Sequence read count per replicate
                                                 # (proxy for sequencing depth)
    # Set Y to NA when sequencing depth is 0 (sequencing failed for all species)
    for (i in 1:I){
      for (j in 1:J){
        for (k in 1:K){
          for (c in 1:C){
            if (M[i, j, k, c] == 0){
              Y[, i, j, k, c] <- NA
            }
          }
        }
      }
    }
  }


  # Extract dimensions names
  names <- list(species = dimnames(Y)[[1]],
                sites = dimnames(Y)[[2]],
                samples = dimnames(Y)[[3]])
  if (protocol != 'seq_read'){
    names <- c(names,
               list(replicates = dimnames(Y)[[4]],
               campaigns = dimnames(Y)[[5]]))
  }else{
    names <- c(names,
               list(campaigns = dimnames(Y)[[4]]))
  }
  names <- c(names,
             list(estimates_sp = c(),
                  estimates = c()))


  # Set the covariates
  if(is.null(covariates)){
    covariates <- new('covarray',
                      psi_cov = NULL,
                      psi_cov_sp = NULL,
                      theta_cov = NULL,
                      theta_cov_sp = NULL,
                      p_cov = NULL,
                      p_cov_sp = NULL,
                      phi_cov = NULL,
                      phi_cov_sp = NULL)
  }
  psi_cov <- covariates@psi_cov
  nb_psi_cov <- sum(tail(dim(psi_cov), 1))           # Number of 'psi_cov' covariates

  psi_cov_sp <- covariates@psi_cov_sp
  nb_psi_cov_sp <- sum(tail(dim(psi_cov_sp), 1))     # Number of 'psi_cov_sp' covariates

  theta_cov <- covariates@theta_cov
  nb_theta_cov <- sum(tail(dim(theta_cov), 1))       # Number of 'theta_cov' covariates

  theta_cov_sp <- covariates@theta_cov_sp
  nb_theta_cov_sp <- sum(tail(dim(theta_cov_sp), 1)) # Number of 'theta_cov_sp' covariates

  p_cov <- covariates@p_cov
  nb_p_cov <- sum(tail(dim(p_cov), 1))               # Number of 'p_cov' covariates

  p_cov_sp <- covariates@p_cov_sp
  nb_p_cov_sp <- sum(tail(dim(p_cov_sp), 1))         # Number of 'p_cov_sp' covariates

  phi_cov <- covariates@phi_cov
  nb_phi_cov <- sum(tail(dim(phi_cov), 1))           # Number of 'phi_cov' covariates

  phi_cov_sp <- covariates@phi_cov_sp
  nb_phi_cov_sp <- sum(tail(dim(phi_cov_sp), 1))     # Number of 'phi_cov_sp' covariates


  # Calculate the number of random effects parameters to record
  # (will be used to iterate across model priors)
  nb_params <-  nb_psi_cov + 1 + nb_theta_cov + 1 # For every protocol estimates for 'psi' and 'theta' covariates will be recorded (slope terms)
                                                  # together with 'alpha_psi' and 'alpha_theta' parameters (intercepts)
  if (protocol == 'PCR_rep'){
    nb_params <- nb_params + 1 + nb_p_cov # Estimates for 'p' covariates will also be recorded
                                          # together with 'alpha_p'
  }
  if (protocol == 'seq_read'){
    nb_params <- nb_params + 1 + nb_phi_cov # Estimates for phi covariates will also be recorded
                                            # together with 'alpha_phi'
  }
  if (protocol == 'PCR_rep_seq_read'){
    nb_params <- nb_params + 1 + nb_p_cov + 1 + nb_phi_cov # Estimates for 'p' and 'phi' covariates will also be recorded
                                                           # together with 'alpha_p' and 'alpha_phi'
  }


  # Set the input data list for JAGS
  data <- list(Y = Y,
               N = N,
               I = I,
               J = J,
               C = C,
               nb_params = nb_params,
               nb_phi_cov = nb_phi_cov,
               nb_p_cov = nb_p_cov,
               nb_theta_cov = nb_theta_cov)

  if (protocol != 'seq_read'){
    data <- c(data,
              list(K = K))
    data <- c(data,
              list(alpha_p_count = 1))
  }else{
    data <- c(data,
              list(alpha_p_count = 0))   # 'p' will not be recorded ('alpha_p_count' will be used when iterating across priors)
  }

  if (protocol != 'PCR_rep'){
    data <- c(data,
              list(M = M))
    data <- c(data,
              list(alpha_phi_count = 1))
  }else{
    data <- c(data,
              list(alpha_phi_count = 0)) # 'phi' will not be recorded ('alpha_phi_count' will be used when iterating across priors)
  }


  # Add all species covariates and their number of random effects to the input data list
  # and record their corresponding names
  if (!is.null(phi_cov_sp)){
    data <- c(data,
              list(nb_phi_cov_sp = nb_phi_cov_sp,
                   phi_cov_sp = phi_cov_sp))
    names$estimates_sp <- c(names$estimates_sp,
                            paste0('beta_phi_sp_', dimnames(phi_cov_sp)[[2]]))
  }

  if (!is.null(p_cov_sp)){
    data <- c(data,
              list(nb_p_cov_sp = nb_p_cov_sp,
                   p_cov_sp = p_cov_sp))
    names$estimates_sp <- c(names$estimates_sp,
                            paste0('beta_p_sp_', dimnames(p_cov_sp)[[2]]))
  }

  if (!is.null(theta_cov_sp)){
    data <- c(data,
              list(nb_theta_cov_sp = nb_theta_cov_sp,
                   theta_cov_sp = theta_cov_sp))
    names$estimates_sp <- c(names$estimates_sp,
                            paste0('beta_theta_sp_', dimnames(theta_cov_sp)[[2]]))
  }

  if (!is.null(psi_cov_sp)){
    data <- c(data,
              list(nb_psi_cov_sp = nb_psi_cov_sp,
                   psi_cov_sp = psi_cov_sp))
    names$estimates_sp <- c(names$estimates_sp,
                            paste0('beta_psi_sp_', dimnames(psi_cov_sp)[[2]]))
  }

  # Add all other covariates and their number of random effects to the input data list
  # and record their corresponding names
  if (protocol != 'PCR_rep'){
    names$estimates <- c(names$estimates, 'alpha_phi')
  }

  if (!is.null(phi_cov) && protocol == 'PCR_rep_seq_read'){
    data <- c(data,
              list(phi_cov = phi_cov))
    names$estimates <- c(names$estimates,
                         paste0('beta_phi_', dimnames(phi_cov)[[5]]))
  }

  if (!is.null(phi_cov) && protocol == 'seq_read'){
    data <- c(data,
              list(phi_cov = phi_cov))
    names$estimates <- c(names$estimates,
                         paste0('beta_phi_', dimnames(phi_cov)[[4]]))
  }

  if (protocol != 'seq_read'){
    names$estimates <- c(names$estimates, 'alpha_p')
  }

  if (!is.null(p_cov)){
    data <- c(data,
              list(p_cov = p_cov))
    names$estimates <- c(names$estimates,
                         paste0('beta_p_', dimnames(p_cov)[[5]]))
  }

  names$estimates <- c(names$estimates, 'alpha_theta')
  if (!is.null(theta_cov)){
    data <- c(data,
              list(theta_cov = theta_cov))
    names$estimates <- c(names$estimates,
                         paste0('beta_theta_', dimnames(theta_cov)[[4]]))
  }

  names$estimates <- c(names$estimates, 'alpha_psi')
  if (!is.null(psi_cov)){
    data <- c(data,
              list(psi_cov = psi_cov,
                   nb_psi_cov = nb_psi_cov)) # 'nb_psi_cov' is added in the list only if psi covariates are added in the model because it will not be used otherwise
                                             # (we do not refer to it when iterating across priors)
    names$estimates <- c(names$estimates,
                         paste0('beta_psi_', dimnames(psi_cov)[[3]]))
  }


  # Write the model in BUGS language
  if (protocol == 'PCR_rep'){
    model <- "
model
{
  # Replicate-level presence/absence (= Detection)
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (k in 1:K){
          for (c in 1:C){
            Y[n, i, j, k, c] ~ dbern(A[n, i, j, c] * p[n, i, j, k, c])"

    if(loglik){
      model <- paste(model, "
            loglik[n, i, j, k, c] <- logdensity.bern(Y[n, i, j, k, c], A[n, i, j, c] * p[n, i, j, k, c])")
    }

    model <- paste(model, "
            logit(p[n, i, j, k, c]) <- alpha_p[n]")

    if (!is.null(p_cov_sp)){
      for(q in 1:nb_p_cov_sp){
        model <- paste(model, paste0("
                                     + beta_p_sp[", q, "] * p_cov_sp[n, ", q, "]"))
      }
    }

    if (!is.null(p_cov)){
      for(q in 1:nb_p_cov){
        model <- paste(model, paste0("
                                     + beta_p[n, ", q, "] * p_cov[i, j, k, c, ", q, "]"))
      }
    }

    model <- paste(model, "
          }
        }
      }
    }
  }")
  }

  if (protocol == 'seq_read'){
    model <- "
model
{
  # Replicate-level sequence read count
  for (i in 1:I){
    for (j in 1:J){
      for (c in 1:C){
        Y[1:N, i, j, c] ~ dmulti(pi[1:N, i, j, c], M[i, j, c])"

    if (loglik){
      model <- paste(model, "
        loglik[i, j, c] <- logdensity.multi(Y[1:N, i, j, c], pi[1:N, i, j, c], M[i, j, c])")
    }

    model <- paste(model, "
      }
    }
  }
  # Sample-level species relative sequence read count
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (c in 1:C){
          ZI[n, i, j, c] ~ dbern(epsilon[n, i, j, c])
          pi[n, i, j, c] <- A[n, i, j, c] * S[n, i, j, c] / ifelse(sum(A[1:N, i, j, c] * S[1:N, i, j, c]) > 0,
                                                                   sum(A[1:N, i, j, c] * S[1:N, i, j, c]),
                                                                   1)
          S[n, i, j, c] ~ dnegbin(delta / (delta + phi[n, i, j, c] * (1 - ZI[n, i, j, c])),
                                  delta)
          log(phi[n, i, j, c]) <- alpha_phi[n]")

    if (!is.null(phi_cov_sp)){
      for(q in 1:nb_phi_cov_sp){
        model <- paste(model, paste0("
                                + beta_phi_sp[", q, "] * phi_cov_sp[n, ", q, "]"))
      }
    }

    if (!is.null(phi_cov)){
      for(q in 1:nb_phi_cov){
        model <- paste(model, paste0("
                                + beta_phi[n, ", q, "] * phi_cov[i, j, c, ", q, "]"))
      }
    }

    model <- paste(model, "
        }
      }
    }
  }")
  }

  if (protocol == 'PCR_rep_seq_read'){
    model <- "
model
{
  # Replicate-level sequence read count
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K) {
        for (c in 1:C) {
          Y[1:N, i, j, k, c] ~ dmulti(pi[1:N, i, j, k, c], M[i, j, k, c])"

    if (loglik){
      model <- paste(model, "
          loglik[i, j, k, c] <- logdensity.multi(Y[1:N, i, j, k, c], pi[1:N, i, j, k, c], M[i, j, k, c])")
      }

    model <- paste(model, "
        }
      }
    }
  }
  # Replicate-level species relative sequence read count
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (k in 1:K){
          for (c in 1:C){
            ZI[n, i, j, k, c] ~ dbern(epsilon[n, i, j, k, c])
            pi[n, i, j, k, c] <- W[n, i, j, k, c] * S[n, i, j, k, c] / ifelse(sum(W[1:N, i, j, k, c] * S[1:N, i, j, k, c]) > 0,
                                                                              sum(W[1:N, i, j, k, c] * S[1:N, i, j, k, c]),
                                                                              1)
            S[n, i, j, k, c] ~ dnegbin(delta / (delta + phi[n, i, j, k, c] * (1 - ZI[n, i, j, k, c])),
                                       delta)
            log(phi[n, i, j, k, c]) <- alpha_phi[n]")

    if (!is.null(phi_cov_sp)){
      for(q in 1:nb_phi_cov_sp){
        model <- paste(model, paste0("
                                     + beta_phi_sp[", q, "] * phi_cov_sp[n, ", q, "]"))
      }
    }

    if (!is.null(phi_cov)){
      for(q in 1:nb_phi_cov){
        model <- paste(model, paste0("
                                     + beta_phi[n, ", q, "] * phi_cov[i, j, k, c, ", q, "]"))
      }
    }

    model <- paste(model, "
          }
        }
      }
    }
  }
  # Replicate-level presence/absence (= Detection)
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (k in 1:K){
          for (c in 1:C){
            W[n, i, j, k, c] ~ dbern(A[n, i, j, c] * p[n, i, j, k, c])
            logit(p[n, i, j, k, c]) <- alpha_p[n]")

    if (!is.null(p_cov_sp)){
      for(q in 1:nb_p_cov_sp){
        model <- paste(model, paste0("
                                     + beta_p_sp[", q, "] * p_cov_sp[n, ", q, "]"))
      }
    }

    if (!is.null(p_cov)){
      for(q in 1:nb_p_cov){
        model <- paste(model, paste0("
                                     + beta_p[n, ", q, "] * p_cov[i, j, k, c, ", q, "]"))
      }
    }

    model <- paste(model, "
          }
        }
      }
    }
  }")
  }

  model <- paste(model, "
  # Sample-level presence/absence (= Collection)
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (c in 1:C){
          A[n, i, j, c] ~ dbern(Z[n, i, c] * theta[n, i, j, c])
          logit(theta[n, i, j, c]) <- alpha_theta[n]")

  if (!is.null(theta_cov_sp)){
    for(q in 1:nb_theta_cov_sp){
      model <- paste(model, paste0("
                                    + beta_theta_sp[", q, "] * theta_cov_sp[n, ", q, "]"))
    }
  }

  if (!is.null(theta_cov)){
    for(q in 1:nb_theta_cov){
      model <- paste(model, paste0("
                                    + beta_theta[n, ", q, "] * theta_cov[i, j, c, ", q, "]"))
    }
  }

  model <- paste(model, "
        }
      }
    }
  }")

  model <- paste(model, "
  # Site-level presence/absence (= Occurrence)
  for (n in 1:N){
    for (i in 1:I){
      for (c in 1:C){
        Z[n, i, c] ~ dbern(psi[n, i, c])
        logit(psi[n, i, c]) <- alpha_psi[n]")

  if (!is.null(psi_cov_sp)){
    for(q in 1:nb_psi_cov_sp){
      model <- paste(model, paste0("
                             + beta_psi_sp[", q, "] * psi_cov_sp[n, ", q, "]"))
    }
  }

  if (!is.null(psi_cov)){
    for(q in 1:nb_psi_cov){
      model <- paste(model, paste0("
                             + beta_psi[n, ", q, "] * psi_cov[i, c, ", q, "]"))
    }
  }

  model <- paste(model, "
      }
    }
  }
")

  model <- paste(model, "
  # Priors")

  if (protocol != 'PCR_rep'){
    model <- paste0(model, "
  delta ~ dgamma(", rho,", ", lambda,")")
  }

  if (!is.null(phi_cov_sp)){
    model <- paste0(model, "
  for (q in 1:nb_phi_cov_sp){
    beta_phi_sp[q] ~ dnorm(0, ", tau,")
  }")
  }

  if (!is.null(p_cov_sp)){
    model <- paste0(model, "
  for (q in 1:nb_p_cov_sp){
    beta_p_sp[q] ~ dnorm(0, ", tau,")
  }")
  }

  if (!is.null(theta_cov_sp)){
    model <- paste0(model, "
  for (q in 1:nb_theta_cov_sp){
    beta_theta_sp[q] ~ dnorm(0, ", tau,")
  }")
  }

  if (!is.null(psi_cov_sp)){
    model <- paste0(model, "
  for (q in 1:nb_psi_cov_sp){
    beta_psi_sp[q] ~ dnorm(0, ", tau,")
  }")
  }

  if (protocol == 'PCR_rep_seq_read'){
    model <- paste0(model, "
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (k in 1:K){
          for (c in 1:C){
            epsilon[n, i, j, k, c] ~ dbeta(1, 1)
          }
        }
      }
    }
  }")
  }

  if (protocol == 'seq_read'){
    model <- paste0(model, "
  for (n in 1:N){
    for (i in 1:I){
      for (j in 1:J){
        for (c in 1:C){
          epsilon[n, i, j, c] ~ dbeta(1, 1)
        }
      }
    }
  }")
  }

  model <- paste(model, "
  for (n in 1:N){")

  if (protocol != 'PCR_rep'){
    model <- paste0(model, "
    alpha_phi[n] <- effect_vector[n, alpha_phi_count]")
  }

  if (!is.null(phi_cov)){
    model <- paste(model, "
    for (q in 1:nb_phi_cov){
      beta_phi[n, q] <- effect_vector[n, alpha_phi_count + q]
    }")
  }

  if (protocol != 'seq_read'){
    model <- paste(model, "
    alpha_p[n] <- effect_vector[n, alpha_phi_count + nb_phi_cov + alpha_p_count]")
  }

  if (!is.null(p_cov)){
    model <- paste(model, "
    for (q in 1:nb_p_cov){
      beta_p[n, q] <- effect_vector[n, alpha_phi_count + nb_phi_cov + alpha_p_count + q]
    }")
  }

    model <- paste(model, "
    alpha_theta[n] <- effect_vector[n, alpha_phi_count + nb_phi_cov + alpha_p_count + nb_p_cov + 1]")

  if (!is.null(theta_cov)){
    model <- paste(model, "
    for (q in 1:nb_theta_cov){
      beta_theta[n, q] <- effect_vector[n, alpha_phi_count + nb_phi_cov + alpha_p_count + nb_p_cov + 1 + q]
    }")
  }

  model <- paste(model, "
    alpha_psi[n] <- effect_vector[n, alpha_phi_count + nb_phi_cov + alpha_p_count + nb_p_cov + 1 + nb_theta_cov + 1]")

  if (!is.null(psi_cov)){
    model <- paste(model, "
    for (q in 1:nb_psi_cov){
      beta_psi[n, q] <- effect_vector[n, alpha_phi_count + nb_phi_cov + alpha_p_count + nb_p_cov + 1 + nb_theta_cov + 1 + q]
    }")
  }

    model <- paste(model, "
    effect_vector[n, 1:nb_params] ~ dmnorm(mu[1:nb_params], omega[1:nb_params, 1:nb_params])
  }
")


    model <- paste(model, "
  # Hyperpriors
  for (q in 1:nb_params){")

    model <- paste0(model, "
    mu[q] ~ dnorm(0, ", tau,")
  }
  omega[1:nb_params, 1:nb_params] ~ dwish(R[1:nb_params, 1:nb_params], nb_params)
  for (q1 in 1:nb_params){
    for (q2 in 1:nb_params){
      R[q1, q2] <- ifelse(q1 == q2, 1, 0)
    }
  }
}")

  writeLines(model, paste0(name, ".txt"))


  # Set the number of lists to initialize
  if (parallel){
    nb_init_lists <- 1 # One list if chains are parallelized
  }else{
    nb_init_lists <- nb_chains # As many lists as chains if they are not parallelized
  }


  # Initialize the model for JAGS
  inits <- list()
  for (c in 1:nb_init_lists){
    init_chain <- list(Z = array(1, dim = c(N, I, C)),
                       A = array(1, dim = c(N, I, J, C)))
    if (protocol == 'PCR_rep_seq_read'){
      init_chain <- c(init_chain,
                      list(W = array(1, dim = c(N, I, J, K, C))),
                      list(S = array(rnbinom(N * I * J * K * C,
                                             size = size,
                                             prob = prob) + 1, # We add 1 to avoid negative values
                                     dim = c(N, I, J, K, C))))
    }
    if (protocol == 'seq_read'){
      init_chain <- c(init_chain,
                      list(S = array(rnbinom(N * I * J * C,
                                             size = size,
                                             prob = prob) + 1, # We add 1 to avoid negative values
                                     dim = c(N, I, J, C))))
    }
    inits[[c]] <- init_chain
  }


  # Set the parameters to record
  parameters <- c('mu',
                  'omega',
                  'alpha_psi',
                  'alpha_theta')

  if (protocol != 'seq_read'){
    parameters <- c(parameters,
                    'alpha_p')
  }
  if (protocol != 'PCR_rep'){
    parameters <- c(parameters,
                    'alpha_phi')
    parameters <- c(parameters,
                    'delta')
  }
  if (!is.null(phi_cov_sp)){
    parameters <- c(parameters,
                    'beta_phi_sp')
  }
  if (!is.null(phi_cov)){
    parameters <- c(parameters,
                    'beta_phi')
  }
  if (!is.null(p_cov_sp)){
    parameters <- c(parameters,
                    'beta_p_sp')
  }
  if (!is.null(p_cov)){
    parameters <- c(parameters,
                    'beta_p')
  }
  if (!is.null(theta_cov_sp)){
    parameters <- c(parameters,
                    'beta_theta_sp')
  }
  if (!is.null(theta_cov)){
    parameters <- c(parameters,
                    'beta_theta')
  }
  if (!is.null(psi_cov_sp)){
    parameters <- c(parameters,
                    'beta_psi_sp')
  }
  if (!is.null(psi_cov)){
    parameters <- c(parameters,
                    'beta_psi')
  }
  if (loglik){
    parameters <- c(parameters,
                    'loglik')
  }
  if ('Z' %in% latent){
    parameters <- c(parameters,
                    'Z')
  }
  if ('psi' %in% posterior){
    parameters <- c(parameters,
                    'psi')
  }
  if ('A' %in% latent){
    parameters <- c(parameters,
                    'A')
  }
  if ('theta' %in% posterior){
    parameters <- c(parameters,
                    'theta')
  }
  if (protocol == 'PCR_rep_seq_read'){
    if ('W' %in% latent){
      parameters <- c(parameters,
                      'W')
    }
  }
  if (protocol != 'seq_read'){
    if ('p' %in% posterior){
      parameters <- c(parameters,
                      'p')
    }
  }
  if (protocol != 'PCR_rep'){
    if ('S' %in% latent){
      parameters <- c(parameters,
                      'S')
    }
    if ('pi' %in% posterior){
      parameters <- c(parameters,
                      'pi')
    }
    if ('phi' %in% posterior){
      parameters <- c(parameters,
                      'phi')
    }
  if ('Y' %in% latent){
    parameters <- c(parameters,
                    'Y')
    }
  }


  # Run JAGS
  if (parallel){
    occupancy_model <- jags.parallel(
      data = data,
      inits = inits,
      parameters.to.save = parameters,
      model.file = paste0(name, ".txt"),
      n.burnin = nb_burnin,
      n.iter = nb_iterations,
      n.thin = nb_thinning,
      n.chains = nb_chains)
  }else{
    occupancy_model <- jags(
      data = data,
      inits = inits,
      parameters.to.save = parameters,
      model.file = paste0(name, ".txt"),
      n.burnin = nb_burnin,
      n.iter = nb_iterations,
      n.thin = nb_thinning,
      n.chains = nb_chains)
  }


  # Return a class object
  instance <- new('Nemodel',
                  model = occupancy_model,
                  protocol = protocol,
                  array = array,
                  covariates = covariates,
                  names = names,
                  loglik = loglik)

  return(instance)
}
