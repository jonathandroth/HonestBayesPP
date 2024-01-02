#' @description
#' This function calculates estimates of the parameters mu and sigma^2 assuming 
#' that the violation of parallel trends between periods, w = delta_t - delta_{t-1}
#' is iid N(mu,sigma^2), i.e. that w_t follows a random walk with normal innovations
#' Given pre-treatment estimates, we compute the maximum likelihood estimates of w and sigma^2

eb_params_ar1 <- function(w, SigmaW){
  
  k <- length(w)
  onevec <- matrix(1, nrow = k, ncol =1)
  A_sigmasq <- function(sigmasq){ return( SigmaW + sigmasq * diag(k)) }
  Ainv_sigmasq <- function(sigmasq){ return(solve(A_sigmasq(sigmasq)) ) }
  
  mu_sigmasq <- function(sigmasq){ return(as.numeric(t(onevec) %*% Ainv_sigmasq(sigmasq) %*% w / as.numeric(t(onevec) %*% Ainv_sigmasq(sigmasq) %*% onevec ) )) }
  
  foc_sigmasq <- function(sigmasq){
    A <- A_sigmasq(sigmasq)
    Ainv <- Ainv_sigmasq(sigmasq)
    mu <- mu_sigmasq(sigmasq)
    muvec <- mu * onevec
    
    foc <-
    t(onevec) %*% ( (A - (w - muvec) %*% t(w - muvec) ) * (-Ainv %*% Ainv )   ) %*% onevec
    
    return(as.numeric(foc))
  }
  
  if(foc_sigmasq(0) < 0){
    sigmasq <- 0
  }else{
    naive_sigmasq <- var(w)
    sigmasq <- stats::uniroot(f=foc_sigmasq, interval = c(0, 10^4*naive_sigmasq))$root
  }
  
  mu <- mu_sigmasq(sigmasq)
  
  return(list(mu = mu, sigmasq = sigmasq))
}


#XX this assumes that betaPre is ordered and the referencePeriod is afterwards
calc_w <- 
function(beta, Sigma, prePeriodIndices){
  
  betaPre <- beta[prePeriodIndices]
  sigmaPre <- Sigma[prePeriodIndices, prePeriodIndices]
  
  numPrePeriods <- length(prePeriodIndices)
  
  #Create the matrix Wmat such that Wmat %*% betapre = w
  # Where w is vector of first-difs
  Wmat <- matrix(data = 0, 
                 nrow = numPrePeriods,
                 ncol = numPrePeriods)
  
  for(t in 1:numPrePeriods){
    #Weight of -1 on current period
    Wmat[t,t] <- -1
    
    #If we're not in the latest period, put weight 1 on next period
    # This gives us delta_t - delta_{t-1} for all but the end
    # For -1, this gives us just -delta_{-1}, since delta_0 = 0
    if(t != numPrePeriods){
      Wmat[t,t+1] <- 1
    }
  }
  
  w <- Wmat %*% betaPre
  SigmaW <- Wmat %*% sigmaPre %*% t(Wmat)
  
  return(list(w = w,
              SigmaW = SigmaW))
}
