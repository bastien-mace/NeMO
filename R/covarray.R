#' @include classes.R
#'
#' @title Format Covariates for NeMO
#'
#' @description This function facilitates the integration of covariates into the
#' **NeMO** modelling framework by formatting them into arrays compatible with [Nemodel].
#'
#' @param protocol Character string. Specifies the modelling protocol to be used in [Nemodel]. \cr
#'   Options are:
#'   - `'PCR_rep'`: Requires a 5D (N × I × J × K × C) presence/absence array
#'   - `'seq_read'`: Requires a 4D (N × I × J × C) sequence read count array
#'   - `'PCR_rep_seq_read'`: Requires a 5D (N × I × J × K × C) sequence read count array \cr
#'   where:
#'      - N: Number of species
#'      - I: Number of sites
#'      - J: Number of samples
#'      - K: Number of PCR replicates
#'      - C: Number of campaigns
#' @param array Input data array to be used in [Nemodel]. Its required dimensionality depends on the selected `protocol` argument.
#' @param cov_list A list of sublists. Each sublist represent a covariate. Sublists must be attributed a unique name and contain the following components:
#'   - `cov_data`: Array of covariate values. These values should be either Boolean for categorical/semi-quantitative predictors or standardised values for quantitative predictors.
#'   - `level`: Character string. Specifies the hierarchical level where the covariate applies (*e.g.*, `'psi'`, `'theta'`, `'p'`, `'phi'`).
#'   - `dimension`: Character string. Specifies the dimensions associated to the covariate. Acceptable values include:
#'       - Single dimensions: `'species'`, `'site'`, `'sample'`, `'replicate'`, or `'campaign'`.
#'       - Combined dimensions: Combinations of the above, except `'species'` (always stands alone). In combinations, terms must be separated by an underscore `'_'`, following this order: 1) `'site'` 2) `'sample'` 3) `'replicate'` 4) `'campaign'` (*i.e.*, `'site_campaign'`).
#'
#' @details The `cov_data` array in each sublist must align with the specified `dimension`. For instance, if `dimension = 'site_sample_campaign'`, the covariate array should be structured such that dimension 1 corresponds to sites, dimension 2 corresponds to samples, and dimension 3 corresponds to campaigns.
#'
#' @return A structured class object with eight slots, organizing covariates across hierarchical levels for direct implementation in the NeMO modelling framework:
#'   - `psi_cov`: Covariates for spatial/methodological/temporal components of \eqn{\psi}.
#'   - `psi_cov_sp`: Species covariates of \eqn{\psi}.
#'   - `theta_cov`: Covariates for spatial/methodological/temporal components of \eqn{\theta}.
#'   - `theta_cov_sp`: Species covariates of \eqn{\theta}.
#'   - `p_cov`: Covariates for spatial/methodological/temporal components of \eqn{p}.
#'   - `p_cov_sp`: Species covariates of \eqn{p}.
#'   - `phi_cov`: Covariates for spatial/methodological/temporal components of \eqn{\varphi}.
#'   - `phi_cov_sp`: Species covariates of \eqn{\varphi}.
#'
#' @examples
#' # Load input array
#' data(fish_PCR_rep)
#'
#' # Load covariate data (Distance to sea)
#' data(distance_cov)
#'
#' # Build the covariate array applied to sites on the psi level
#' covarray(protocol = 'PCR_rep',
#'          array = fish_PCR_rep,
#'          cov_list = list(Distance = list(cov_data = distance_cov$Distance,
#'                                          level = 'psi',
#'                                          dimension = 'site')))
#'
#' @seealso \code{\link{Nemodel}}
#'
#' @export
#'
#' @importFrom abind abind
#' @importFrom magrittr %>%
#' @importFrom methods new
#'
covarray <- function(protocol = c('PCR_rep', 'seq_read', 'PCR_rep_seq_read'),
                     array,
                     cov_list = list(
                        list(cov_data = NULL,
                             level = 'psi',
                             dimension = 'species')))
{

  # Ensure protocol is provided and valid
  if (missing(protocol)){
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

  # Define the allowed values for 'level' and 'dimension'
  allowed_levels <- c('psi', 'theta', 'p', 'phi')
  allowed_dimensions <- c('species', 'site', 'sample', 'replicate', 'campaign',
                          'site_sample', 'site_replicate', 'site_campaign', 'sample_replicate', 'sample_campaign', 'replicate_campaign',
                          'site_sample_replicate', 'site_sample_campaign', 'site_replicate_campaign', 'sample_replicate_campaign',
                          'site_sample_replicate_campaign')

  # Check 'cov_list' components
  for(name in names(cov_list)){
    if(is.null(name) || name == ''){
      stop("Error: Please provide a name to your lists.")
      }
    }
  for (cov in cov_list){
    if(sum(is.na(cov$cov_data)) > 0){
      stop("Error: NA values are not allowed.")
    }
    if (!is.list(cov) || !all(c('level', 'dimension', 'cov_data') %in% names(cov))){
      stop("Error: Each element of 'cov_list' must be a list with 'level', 'dimension', and 'cov_data'.")
    }
    if (!cov$level %in% allowed_levels){
      stop("Error: 'level' in 'cov_list' can only contain the following values: ",
           paste(allowed_levels, collapse = ', '))
    }
    if (!cov$dimension %in% allowed_dimensions){
      stop("Error: 'dimension' in 'cov_list' can only contain the following values: ",
           paste(allowed_dimensions, collapse = ', '))
    }
  }

  # Check that cov_data is numeric (integer, double)
  for (cov in cov_list) {
    if (!any(sapply(cov$cov_data, is.numeric))) {
      stop("Error: Values in 'cov_data' must be numeric.")
    }
  }


  # Set input array dimensions
  N <- dim(array)[1] # Species
  I <- dim(array)[2] # Sites
  J <- dim(array)[3] # Samples
  if (protocol != 'seq_read'){
    K <- dim(array)[4] # Replicates
    C <- dim(array)[5] # Campaigns
  }else{
    C <- dim(array)[4]  # Campaigns
  }


  # Initialize covariates arrays
  psi_cov <- NULL
  psi_cov_sp <- NULL
  theta_cov <- NULL
  theta_cov_sp <- NULL
  p_cov <- NULL
  p_cov_sp <- NULL
  phi_cov <- NULL
  phi_cov_sp <- NULL


  # Count the number of covariates associated to each levels
  if (!(is.null(cov_list))){
    nb_psi <- sum(sapply(cov_list, function(x) x$level == 'psi'))
    nb_theta <- sum(sapply(cov_list, function(x) x$level == 'theta'))
    nb_p <- sum(sapply(cov_list, function(x) x$level == 'p'))
    nb_phi <- sum(sapply(cov_list, function(x) x$level == 'phi'))
  }else{
    nb_psi <- 0
    nb_theta <- 0
    nb_p <- 0
    nb_phi <- 0
  }


  # Build the covariate arrays associated with 'psi' level
  if (nb_psi > 0){
    for (i in 1:length(cov_list)){
      cov <- cov_list[[i]]

      # Build the 'psi_cov_sp' covariate array (Species covariate)
      if (cov$level == 'psi' && cov$dimension == 'species'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != N){
          stop("Error: 'cov_data' must be of length ", N, " when 'dimension' is species.")
        }
        psi_sp_array <- array(unlist(cov$cov_data[]), dim = c(N)) # Give correct dimensions
        psi_cov_sp <- abind(psi_cov_sp, psi_sp_array, along = 2) # Add covariate values to a new layer
        dimnames(psi_cov_sp)[[2]][dim(psi_cov_sp)[[2]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }

      # Build the 'psi_cov' covariate array
      if (cov$level == 'psi' && cov$dimension == 'site'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != I){
          stop("Error: 'cov_data' must be of length ", I, " when 'dimension' is site.")
        }
        psi_array <- array(unlist(rep(cov$cov_data[], C)), dim = c(I, C)) # Give correct dimensions
        psi_cov <- abind(psi_cov, psi_array, along = 3) # Add covariate values to a new layer
        dimnames(psi_cov)[[3]][dim(psi_cov)[[3]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'psi' && cov$dimension == 'campaign'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != C){
          stop("Error: 'cov_data' must be of length ", C, " when 'dimension' is campaign.")
        }
        psi_array <- array(unlist(rep(cov$cov_data[], I)), dim = c(C, I)) %>%
          aperm(c(2, 1)) # Give correct dimensions
        psi_cov <- abind(psi_cov, psi_array, along = 3) # Add covariate values to a new layer
        dimnames(psi_cov)[[3]][dim(psi_cov)[[3]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'psi' && cov$dimension == 'site_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * C, " when 'dimension' is site_campaign.")
        }
        psi_array <- array(unlist(rep(cov$cov_data[,])), dim = c(I, C)) # Give correct dimensions
        psi_cov <- abind(psi_cov, psi_array, along = 3) # Add covariate values to a new layer
        dimnames(psi_cov)[[3]][dim(psi_cov)[[3]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
    }
  }


  # Build the covariate arrays associated with 'theta' level
  if (nb_theta > 0){
    for (i in 1:length(cov_list)){
      cov <- cov_list[[i]]

      # Build the 'theta_cov_sp' covariate array (Species covariate)
      if (cov$level == 'theta' && cov$dimension == 'species'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != N){
          stop("Error: 'cov_data' must be of length ", N, " when 'dimension' is species.")
        }
        theta_sp_array <- array(unlist(cov$cov_data[]), dim = c(N)) # Give correct dimensions
        theta_cov_sp <- abind(theta_cov_sp, theta_sp_array, along = 2) # Add covariate values to a new layer
        dimnames(theta_cov_sp)[[2]][dim(theta_cov_sp)[[2]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }

      # Build the 'theta_cov' covariate array
      if (cov$level == 'theta' && cov$dimension == 'site'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != I){
          stop("Error: 'cov_data' must be of length ", I, " when 'dimension' is site.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[], J * C)), dim = c(I, J, C)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'theta' && cov$dimension == 'sample'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != J){
          stop("Error: 'cov_data' must be of length ", J, " when 'dimension' is sample.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[], I * C)), dim = c(J, I, C)) %>%
          aperm(c(2, 1, 3)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'theta' && cov$dimension == 'campaign'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != C){
          stop("Error: 'cov_data' must be of length ", C, " when 'dimension' is campaign.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[], I * J)), dim = c(C, I, J)) %>%
          aperm(c(2, 3, 1)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'theta' && cov$dimension == 'site_sample_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * J * C, " when 'dimension' is site_sample_campaign.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[,,])), dim = c(I, J, C)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'theta' && cov$dimension == 'site_sample'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J){
          stop("Error: 'cov_data' must be of dimensions ", I * J, " when 'dimension' is site_sample.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[,], C)), dim = c(I, J, C)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'theta' && cov$dimension == 'site_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * C, " when 'dimension' is site_campaign.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[,], J)), dim = c(I, C, J)) %>%
          aperm(c(1, 3, 2)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'theta' && cov$dimension == 'sample_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != C){
          stop("Error: 'cov_data' must be of dimensions ", J * C, " when 'dimension' is sample_campaign.")
        }
        theta_array <- array(unlist(rep(cov$cov_data[,], I)), dim = c(J, C, I)) %>%
          aperm(c(3, 1, 2)) # Give correct dimensions
        theta_cov <- abind(theta_cov, theta_array, along = 4) # Add covariate values to a new layer
        dimnames(theta_cov)[[4]][dim(theta_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
    }
  }


  # Build the covariate arrays associated with 'p' level
  if (nb_p > 0){
    for (i in 1:length(cov_list)){
      cov <- cov_list[[i]]

      # Build the 'p_cov_sp' covariate array (Species covariate)
      if (cov$level == 'p' && cov$dimension == 'species'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != N){
          stop("Error: 'cov_data' must be of length ", N, " whenn 'level' is psi and 'dimension' is species.")
        }
        p_sp_array <- array(unlist(cov$cov_data[]), dim = c(N)) # Give correct dimensions
        p_cov_sp <- abind(p_cov_sp, p_sp_array, along = 2) # Add covariate values to a new layer
        dimnames(p_cov_sp)[[2]][dim(p_cov_sp)[[2]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }

      # Build the 'p_cov' covariate array
      if (cov$level == 'p' && cov$dimension == 'site'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != I){
          stop("Error: 'cov_data' must be of length ", I, " when 'dimension' is site.")
        }
        p_array <- array(unlist(rep(cov$cov_data[], J * K * C)), dim = c(I, J, K, C)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'sample'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != J){
          stop("Error: 'cov_data' must be of length ", J, " when 'dimension' is sample.")
        }
        p_array <- array(unlist(rep(cov$cov_data[], I * K * C)), dim = c(J, I, K, C)) %>%
          aperm(c(2, 1, 3, 4)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'replicate'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != K){
          stop("Error: 'cov_data' must be of length ", K, " when 'dimension' is replicate.")
        }
        p_array <- array(unlist(rep(cov$cov_data[], I * J * C)), dim = c(K, I, J, C)) %>%
          aperm(c(2, 3, 4, 1)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'campaign'){
        # Check if covariate values were provided with correct dimensions
        if(length(cov$cov_data) != C){
          stop("Error: 'cov_data' must be of length ", C, " when 'dimension' is campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[], I * J * K)), dim = c(C, I, J, K)) %>%
          aperm(c(2, 3, 4, 1)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_sample_replicate_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != K && dim(cov$cov_data)[4] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * J * K * C, " when 'dimension' is site_sample_replicate_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,,,])), dim = c(I, J, K, C)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_sample_replicate'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != K){
          stop("Error: 'cov_data' must be of dimensions ", I * J * K, " when 'dimension' is site_sample_replicate.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,,], C)), dim = c(I, J, K, C)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_sample_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * J * C, " when 'dimension' is site_sample_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,,], N * K)), dim = c(I, J, C, K)) %>%
          aperm(c(1, 2, 4, 3)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_replicate_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != K && dim(cov$cov_data)[3] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * K * C, " when 'dimension' is site_replicate_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,,], J)), dim = c(I, K, C, J)) %>%
          aperm(c(1, 4, 2, 3)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'sample_replicate_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != K && dim(cov$cov_data)[3] != C){
          stop("Error: 'cov_data' must be of dimensions ", J * K * C, " when 'dimension' is sample_replicate_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,,], I)), dim = c(J, K, C, I)) %>%
          aperm(c(4, 5, 1, 2, 3)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_sample'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J){
          stop("Error: 'cov_data' must be of dimensions ", I * J, " when 'dimension' is site_sample.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,], K * C)), dim = c(I, J, K, C)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_replicate'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != K){
          stop("Error: 'cov_data' must be of dimensions ", I * K, " when 'dimension' is site_replicate.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,], J * C)), dim = c(I, K, J, C)) %>%
          aperm(c(1, 3, 2, 4)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'site_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != C){
          stop("Error: 'cov_data' must be of dimensions ", I * C, " when 'dimension' is site_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,], J * K)), dim = c(I, C, J, K)) %>%
          aperm(c(1, 3, 4, 2)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'sample_replicate'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != K){
          stop("Error: 'cov_data' must be of dimensions ", J * K, " when 'dimension' is sample_replicate.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,], I * C)), dim = c(J, K, I, C)) %>%
          aperm(c(3, 1, 2, 4)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'sample_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != C){
          stop("Error: 'cov_data' must be of dimensions ", J * C, " when 'dimension' is sample_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,], I * K)), dim = c(J, C, I, K)) %>%
          aperm(c(3, 1, 4, 2)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
      if (cov$level == 'p' && cov$dimension == 'replicate_campaign'){
        # Check if covariate values were provided with correct dimensions
        if(dim(cov$cov_data)[1] != K && dim(cov$cov_data)[2] != C){
          stop("Error: 'cov_data' must be of dimensions ", K * C, " when 'dimension' is replicate_campaign.")
        }
        p_array <- array(unlist(rep(cov$cov_data[,], I * J)), dim = c(K, C, I, J)) %>%
          aperm(c(3, 4, 1, 2)) # Give correct dimensions
        p_cov <- abind(p_cov, p_array, along = 5) # Add covariate values to a new layer
        dimnames(p_cov)[[5]][dim(p_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
      }
    }
  }


  # Build the covariate arrays associated with 'phi' level (when 'seq_read')
  if (protocol == 'seq_read'){
    if (nb_phi > 0){
      for (i in 1:length(cov_list)){
        cov <- cov_list[[i]]

        # Build the 'phi_cov_sp' covariate array (Species covariate)
        if (cov$level == 'phi' && cov$dimension == 'species'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != N){
            stop("Error: 'cov_data' must be of length ", N, " when 'dimension' is species.")
          }
          phi_sp_array <- array(unlist(cov$cov_data[]), dim = c(N)) # Give correct dimensions
          phi_cov_sp <- abind(phi_cov_sp, phi_sp_array, along = 2) # Add covariate values to a new layer
          dimnames(phi_cov_sp)[[2]][dim(phi_cov_sp)[[2]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }

        # Build the 'phi_cov' covariate array
        if (cov$level == 'phi' && cov$dimension == 'site'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != I){
            stop("Error: 'cov_data' must be of length ", I, " when 'dimension' is site.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], J * C)), dim = c(I, J, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'sample'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != J){
            stop("Error: 'cov_data' must be of length ", J, " when 'dimension' is sample.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], I * C)), dim = c(J, I, C)) %>%
            aperm(c(2, 1, 3)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'campaign'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != C){
            stop("Error: 'cov_data' must be of length ", C, " when 'dimension' is campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], I * J)), dim = c(C, I, J)) %>%
            aperm(c(2, 3, 1)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_sample_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != C){
            stop("Error: 'cov_data' must be of dimensions ", I * J * C, " when 'dimension' is site_sample_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,,])), dim = c(I, J, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_sample'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J){
            stop("Error: 'cov_data' must be of dimensions ", I * J, " when 'dimension' is site_sample.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], C)), dim = c(I, J, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != C){
            stop("Error: 'cov_data' must be of dimensions ", I * C, " when 'dimension' is site_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], J)), dim = c(I, C, J)) %>%
            aperm(c(1, 3, 2)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'sample_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != C){
            stop("Error: 'cov_data' must be of dimensions ", J * C, " when 'dimension' is sample_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], I)), dim = c(J, C, I)) %>%
            aperm(c(3, 1, 2)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 4) # Add covariate values to a new layer
          dimnames(phi_cov)[[4]][dim(phi_cov)[[4]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
      }
    }
  }


  # Build the covariate arrays associated with 'phi' level (when 'PCR_rep_seq_read')
  if (protocol == 'PCR_rep_seq_read'){
    if (nb_phi > 0){
      for (i in 1:length(cov_list)){
        cov <- cov_list[[i]]

        # Build the 'phi_cov_sp' covariate array (Species covariate)
        if (cov$level == 'phi' && cov$dimension == 'species'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != N){
            stop("Error: 'cov_data' must be of length ", N, " whenn 'level' is psi and 'dimension' is species.")
          }
          phi_sphi_array <- array(unlist(cov$cov_data[]), dim = c(N)) # Give correct dimensions
          phi_cov_sp <- abind(phi_cov_sp, phi_sphi_array, along = 2) # Add covariate values to a new layer
          dimnames(phi_cov_sp)[[2]][dim(phi_cov_sp)[[2]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }

        # Build the 'phi_cov' covariate array
        if (cov$level == 'phi' && cov$dimension == 'site'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != I){
            stop("Error: 'cov_data' must be of length ", I, " when 'dimension' is site.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], J * K * C)), dim = c(I, J, K, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'sample'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != J){
            stop("Error: 'cov_data' must be of length ", J, " when 'dimension' is sample.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], I * K * C)), dim = c(J, I, K, C)) %>%
            aperm(c(2, 1, 3, 4)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'replicate'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != K){
            stop("Error: 'cov_data' must be of length ", K, " when 'dimension' is replicate.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], I * J * C)), dim = c(K, I, J, C)) %>%
            aperm(c(2, 3, 4, 1)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'campaign'){
          # Check if covariate values were provided with correct dimensions
          if(length(cov$cov_data) != C){
            stop("Error: 'cov_data' must be of length ", C, " when 'dimension' is campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[], I * J * K)), dim = c(C, I, J, K)) %>%
            aperm(c(2, 3, 4, 1)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_sample_replicate_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != K && dim(cov$cov_data)[4] != C){
            stop("Error: 'cov_data' must be of dimensions ", I * J * K * C, " when 'dimension' is site_sample_replicate_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,,,])), dim = c(I, J, K, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_sample_replicate'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != K){
            stop("Error: 'cov_data' must be of dimensions ", I * J * K, " when 'dimension' is site_sample_replicate.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,,], C)), dim = c(I, J, K, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_sample_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J && dim(cov$cov_data)[3] != C){
            stop("Error: 'cov_data' must be of dimensions ", I * J * C, " when 'dimension' is site_sample_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,,], N * K)), dim = c(I, J, C, K)) %>%
            aperm(c(1, 2, 4, 3)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_replicate_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != K && dim(cov$cov_data)[3] != C){
            stop("Error: 'cov_data' must be of dimensions ", I * K * C, " when 'dimension' is site_replicate_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,,], J)), dim = c(I, K, C, J)) %>%
            aperm(c(1, 4, 2, 3)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'sample_replicate_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != K && dim(cov$cov_data)[3] != C){
            stop("Error: 'cov_data' must be of dimensions ", J * K * C, " when 'dimension' is sample_replicate_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,,], I)), dim = c(J, K, C, I)) %>%
            aperm(c(4, 5, 1, 2, 3)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_sample'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != J){
            stop("Error: 'cov_data' must be of dimensions ", I * J, " when 'dimension' is site_sample.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], K * C)), dim = c(I, J, K, C)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_replicate'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != K){
            stop("Error: 'cov_data' must be of dimensions ", I * K, " when 'dimension' is site_replicate.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], J * C)), dim = c(I, K, J, C)) %>%
            aperm(c(1, 3, 2, 4)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'site_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != I && dim(cov$cov_data)[2] != C){
            stop("Error: 'cov_data' must be of dimensions ", I * C, " when 'dimension' is site_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], J * K)), dim = c(I, C, J, K)) %>%
            aperm(c(1, 3, 4, 2)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'sample_replicate'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != K){
            stop("Error: 'cov_data' must be of dimensions ", J * K, " when 'dimension' is sample_replicate.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], I * C)), dim = c(J, K, I, C)) %>%
            aperm(c(3, 1, 2, 4)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'sample_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != J && dim(cov$cov_data)[2] != C){
            stop("Error: 'cov_data' must be of dimensions ", J * C, " when 'dimension' is sample_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], I * K)), dim = c(J, C, I, K)) %>%
            aperm(c(3, 1, 4, 2)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
        if (cov$level == 'phi' && cov$dimension == 'replicate_campaign'){
          # Check if covariate values were provided with correct dimensions
          if(dim(cov$cov_data)[1] != K && dim(cov$cov_data)[2] != C){
            stop("Error: 'cov_data' must be of dimensions ", K * C, " when 'dimension' is replicate_campaign.")
          }
          phi_array <- array(unlist(rep(cov$cov_data[,], I * J)), dim = c(K, C, I, J)) %>%
            aperm(c(3, 4, 1, 2)) # Give correct dimensions
          phi_cov <- abind(phi_cov, phi_array, along = 5) # Add covariate values to a new layer
          dimnames(phi_cov)[[5]][dim(phi_cov)[[5]]] <- names(cov_list)[i] # Attribute covariate name to the layer
        }
      }
    }
  }


  # Return a class object
  instance <- new('covarray',
                  psi_cov = psi_cov,           # Dimension: I x C
                  psi_cov_sp = psi_cov_sp,     # Dimension: N
                  theta_cov = theta_cov,       # Dimension: I x J x C
                  theta_cov_sp = theta_cov_sp, # Dimension: N
                  p_cov = p_cov,               # Dimension: I x J x K x C
                  p_cov_sp = p_cov_sp,         # Dimension: N
                  phi_cov = phi_cov,           # Dimension: I x J x C ('PCR_rep') ou I x J x K x C ('PCR_rep_seq_read)
                  phi_cov_sp = phi_cov_sp)     # Dimension: N

  return(instance)
}
