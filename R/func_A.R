#' A = (Y_i - mu_i)
#'
#'
#' @param y binary reponse
#' @param mu mean of the binary response
#' @return A

func_A <- function(y, mu){
  return(y - mu)
}

