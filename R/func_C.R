#' C = \partial E(Yij) / \partial beta = x*mu*(1-mu)
#'
#'
#' @param x covariates
#' @param mu mean of the binary response
#' @return C


func_C <- function(x, mu){
  return(x * matrix(rep(mu,ncol(x)),ncol=ncol(x)) * (1 - matrix(rep(mu,ncol(x)),ncol=ncol(x))))
}
