#' nu_{ijk} = E(Y_{ij} Y_{ik})
#' nu: a choose(ni,2) vector of nu_{ijk}
#'
#' @param gamma a vector of log odds ration within cluster
#' @param mu mean of the binary response
#' @return nu

func_nu <- function(gamma, mu){
  # theta_{ijk}
  theta <- exp(gamma)
  nu <- f_nu(theta, mu)
  return(nu)
}