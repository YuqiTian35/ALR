#' initial values for alpha and beta
#'
#'
#' @param y: binary response
#' @param x: covariates
#' @param q: number of alphas (q=1 for exchangeable)
#' @return A list of initial values for beta and alpha
#' @export 
initial_val <- function(mean.formula, data, q=1){
  # alpha: global odds ratio parameters
  alpha <- rep(0.01, q)
  # beta: mean parameter (initial beta)
  beta <- lrm(mean.formula, data)$coefficients
  return(list(beta=beta, alpha=alpha))
}