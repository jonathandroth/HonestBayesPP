#' @description This function creates the covariance matrix V_delta associated with an AR(1) for the treatment-control dif in each period
#' @param rho The autocorrelation of the AR(1)
#' @param sigma The sd of the AR(1) innovations
#' @param tVec A vector corresponding to the time periods of the event-study coefs
#' @param referencePeriod (Optional) The time period of the omitted coef (default 0)
#' @return A variance-covariance matrix with dimension length(tVec)
create_V_AR1 <- function(rho,sigma,tVec, referencePeriod = 0){
  
  relativeTime <- tVec - referencePeriod
  numPeriods <- length(tVec)
  V <- matrix(nrow = numPeriods,
              ncol = numPeriods)
  
  if(rho != 1){
    #AR(1) case
    for(i in 1:numPeriods){
      for(j in i:numPeriods){
        t <- relativeTime[i]
        tprime <- relativeTime[j]
        
        vij <- (rho^abs(t-tprime) - rho^abs(t) - rho^abs(tprime) +1)/(1-rho^2) * sigma^2
        V[i,j] <- vij
        
        if(i != j){
          V[j,i] <- vij
        }
      }
    }
  }else{
    #Random walk case
    for(i in 1:numPeriods){
      for(j in i:numPeriods){
        t <- relativeTime[i]
        tprime <- relativeTime[j]
        
        if(sign(t) != sign(tprime)){
          vij <- 0
        }else{
          vij <- pmin( abs(t), abs(tprime) ) * sigma^2
        }
        
        V[i,j] <- vij
        if(i != j){
          V[j,i] <- vij
        }
        
      }
    }  
  }  
  return(V)
}

# create_V_AR1(rho = 0.766,
#              sigma = 0.04,
#              tVec = c(seq(-2,-1),seq(1,3)))