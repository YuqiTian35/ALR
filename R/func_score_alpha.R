#' second estimating equation (for alpha)
#'
#'
#' @param beta: initial values of beta
#' @param alpha: initial values of alpha
#' @param y: binary response
#' @param x covariates
#' @param z covariates that speicifies the form of assciation within cluster
#' @param n a vector of sample sizes in each cluster
#' @param q number of pairs/OR parameters (alpha)
#' @return A list of matrices for solving the second estimating equation

func_score_alpha <- function(beta, alpha, y, x, z, n, q){
  # initial
  U <- rep(0, q)
  Ustar <- matrix(0, nrow = q, ncol = q)
  UUtran <- matrix(0, nrow = q, ncol = q)
  
  # start and end points for x and z
  ind_x <- start_end_ind(n)
  ind_z <- start_end_ind(n*(n-1)/2) # n choose 2
  
  for(i in 1:length(n)){ # for each subject
    
    # data in each cluster
    x_c <- x[ind_x$start[i]:ind_x$end[i],]
    y_c <- y[ind_x$start[i]:ind_x$end[i]]
    
    if(n[i] == 1){ # only 1 obs in the cluster
      next
    }else{ # more than 1 obs in the cluster
      if(q == 1){
        z_c <- matrix(z[ind_z$start[i]:ind_z$end[i],], ncol=1)
      }else{
        z_c <- data.frame.to_matrix(z[ind_z$start[i]:ind_z$end[i],])
        if(n[i]==2){
          z_c <- matrix(z_c,nrow=1)
        }
      }
      
      # initialization
      U_c <- rep(0, q)  #score
      Ustar_c <-  matrix(data = 0, nrow = q, ncol = q) # information
      mu_c <- 1 / (1 + exp(-x_c %*% beta))
      gamma_c <- z_c %*% alpha
      
      # elements in 2nd order estimating equation
      second_order <- func_second(y_c, z_c, mu_c, n[i], gamma_c)
      
      t <- second_order$t
      inv_s <- second_order$inv_s
      r <- as.matrix(second_order$r)
      
      # 2nd estimating equation
      U_c <- mat_3(t(t), inv_s, r)
      Ustar_c <- mat_3(t(t), inv_s, t)
      UUtran_c <- mat_sqr(U_c)
      
      U <- U + U_c
      Ustar <- Ustar + Ustar_c
      UUtran <- UUtran + UUtran_c
      
    }
  }
  return(list(U=U, Ustar=Ustar, UUtran=UUtran))
}
