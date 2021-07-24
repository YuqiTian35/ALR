#' first estimating equation (for beta)
#'
#'
#' @param beta: initial values of beta
#' @param alpha: initial values of alpha
#' @param y: binary response
#' @param x covariates
#' @param z covariates that speicifies the form of assciation within cluster
#' @param n a vector of sample sizes in each cluster
#' @param p number of mean parameters (beta)
#' @param q number of pairs/OR parameters (alpha)
#' @return A list of matrices for solving the first estimating equation

func_score_beta <-function(beta, alpha, y, x, z, n, p, q){
  # initial
  U <- rep(0, p)
  Ustar <- matrix(0, nrow = p, ncol = p)
  UUtran <- matrix(0, nrow = p, ncol = p)
  
  # start and end points for x and z
  ind_x <- start_end_ind(n)
  ind_z <- start_end_ind(n*(n-1)/2) # n choose 2
  
  for(i in 1:length(n)){ # for each subject
    
    # data in each cluster
    x_c <- x[ind_x$start[i]:ind_x$end[i],]
    y_c <- y[ind_x$start[i]:ind_x$end[i]]
    
    
    if(n[i] == 1){ # only 1 obs in the cluster
      # x_c <- as.matrix(x_c)
      mu_c <- 1 / (1 + exp(- beta %*% x_c))
      U_c <- x_c %*% (y_c - mu_c) # score
      Ustar_c <- x_c %*% (mu_c * (1 - mu_c)) %*% t(x_c) # information
      
      # update
      U <- U + U_c
      Ustar <- Ustar + Ustar_c
      
    }else{ # more than 1 obs in the cluster
      
      if(q == 1){ # if only 1 OR parameter
        z_c <- matrix(z[ind_z$start[i]:ind_z$end[i],], ncol=1)
      }else{
        z_c <- data.frame.to_matrix(z[ind_z$start[i]:ind_z$end[i],])
        if(n[i]==2){ 
          z_c <- matrix(z_c,nrow=1)
        }
      }
      # initialization
      U_c <- rep(0, p)  #score
      Ustar_c <-  matrix(0, nrow = p, ncol = p) # information
      mu_c <- 1 / (1 + exp(-x_c %*% beta))
      gamma_c <- z_c %*% alpha
      
      # elements in 1st order estimating equation
      c <- func_C(x_c, mu_c)
      nu_c <- func_nu(gamma_c, mu_c)
      b <- func_B(nu_c, mu_c)
      a <- func_A(y_c, mu_c)
      
      solve_b <- mat_inv(b)
      # 1st estimating equation
      U_c <- mat_3(t(c), solve_b, a)
      Ustar_c <- mat_3(t(c), solve_b, c)
      UUtran_c <- mat_sqr(U_c)
      
      U <- U + U_c
      Ustar <- Ustar + Ustar_c
      UUtran <- UUtran + UUtran_c
    }
  }
  return(list(U=U, Ustar=Ustar, UUtran=UUtran))
}