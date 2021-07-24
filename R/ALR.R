#' ALR: Alternating Logistic Regression
#'
#'
#' @param formula a R formula object
#' @param data a data frame
#' @param z a matrix of covariates that specifies the form of assciation between Yij and Yik
#' @param n a vector of sample sizes in each cluster
#' @param maxit number of maximum iteration
#' @param tol tolerace limimt
#' @return A list of results: beta, alpha, Cov(beta), Cov(alpha), number of itereation
#' @expor
ALR <- function(mean.formula, data, z, n, maxit=100, tol=1e-5){

  mean.f <- model.frame(mean.formula, data)

  # data
  y <- model.response(mean.f, "numeric")
  x <- model.matrix(mean.formula, mean.f)

  # p: number of mean parameters
  p <- ncol(x)
  # number of pairs/OR parameters
  q <- ncol(z)

  # update for each parameter
  delta_beta <- rep(2*tol, p)
  delta_alpha <- rep(2*tol, q)

  # initial values for parameters
  i_para <- initial_val(mean.formula, data, q)
  beta <- i_para$beta
  alpha <- i_para$alpha

  # loops
  for(niter in 1:maxit){

    if(max(abs(c(delta_beta, as.numeric(delta_alpha)))) < tol){
      break
    }

    # beta update
    score_beta <- func_score_beta(beta, alpha, y, x, z, n, p, q)
    delta_beta <- solve(score_beta$Ustar, score_beta$U)
    beta <- beta + c(delta_beta)

    # alpha update
    score_alpha <- func_score_alpha(beta, alpha, y, x, z, n, q)
    delta_alpha <- solve(score_alpha$Ustar, score_alpha$U)
    alpha <- alpha + c(delta_alpha)
  }

  # covaraince matrix - beta
  inv_Ustar_beta <- solve(score_beta$Ustar)
  cov_beta <- inv_Ustar_beta %*% score_beta$UUtran %*% inv_Ustar_beta

  # covaraince matrix - alpha
  inv_Ustar_alpha <- solve(score_alpha$Ustar)
  cov_alpha <- inv_Ustar_alpha %*% score_alpha$UUtran %*% inv_Ustar_alpha


  return(list(beta=beta, alpha=alpha,
              cov_beta=cov_beta, cov_alpha = cov_alpha,
              niter=niter))
}
