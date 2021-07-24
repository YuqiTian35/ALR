#' second estmiating equation
#'
#'
#' @param y binary response
#' @param z covariates that speicifies the form of assciation within cluster
#' @param mu mean of the binary response
#' @param n a vector of sample sizes in each cluster
#' @param gamma a vector of log odds ration within cluster
#' @return A list of initial values for beta and alpha


func_second <- function(y, z, mu, n, gamma){
  # l <- 1
  n_choose_2 <- as.integer(choose(n, 2))
  ncol_z <- ncol(z)
  theta <- exp(gamma)
  res <- f_2(y, n, z, mu, gamma, theta, ncol_z)
  return(res)
  
}