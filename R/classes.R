setClass('rjags')
setClassUnion('array_NULL', c('array', 'NULL'))
setClassUnion('list_NULL', c('list', 'NULL'))

setClass('covarray',
         slots = c(
           psi_cov = 'array_NULL',
           psi_cov_sp = 'array_NULL',
           theta_cov = 'array_NULL',
           theta_cov_sp = 'array_NULL',
           p_cov = 'array_NULL',
           p_cov_sp = 'array_NULL',
           phi_cov = 'array_NULL',
           phi_cov_sp = 'array_NULL')
)

setClass('Nemodel',
         slots = c(
           model = 'rjags',
           protocol = 'character',
           array = 'array',
           covariates = 'covarray',
           names = 'list',
           loglik = 'logical')
)

setClass('min_resources',
         slots = c(
           J_min = 'list_NULL',
           K_min = 'list_NULL',
           M_min = 'list_NULL')
)

setClass('WAIC',
         slots = c(
           waic = 'numeric',
           lppd = 'numeric',
           p_waic = 'numeric')
)
