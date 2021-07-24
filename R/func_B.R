#' B = Cov(Y_i) = sqrt(Si) times Ri times sqrt(Si)
#' Cov(Y_ijk) = E(Y_ij Y_ik) - E(Y_ij)E(Y_ik) = nu_{ijk}-mu_{ij}mu_{ik}
#' @param nu E(YijYik)  
#' @param gamma a vector of log odds ration within cluster
#' @param mu mean of the binary response
#' @return B


func_B <- function(nu, mu){
  b <- f_b(nu, mu)
  return(b)
}
