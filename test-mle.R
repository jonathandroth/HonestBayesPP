w <-
calc_w(beta = -3:3,
       Sigma = diag(7),
       prePeriodIndices = 1:3)

expect_equal(w$w,
             matrix(1,nrow=3))

expect_equal(diag(w$SigmaW),
             c(2,2,1))


#Test that MLE works in example with SigmaW diagonal
numPeriods <- 1000
mu <- 4
SigmaW <- diag(numPeriods)
what <- MASS::mvrnorm(n = 1,
                      mu = rep(mu, numPeriods),
                      Sigma = diag(numPeriods) + SigmaW )

mle <- eb_params_ar1(w=what,SigmaW = SigmaW)

expect_equal(mle$mu,mu, tol = 5/sqrt(numPeriods))
expect_equal(mle$sigmasq,1, tol = 5/sqrt(numPeriods))


#Test that MLE works in example with SigmaW non-diagonal SigmaW corresponding to correlation matrix of AR(1)
numPeriods <- 1000
mu <- 4

SigmaW <- diag(numPeriods)
for(i in 1:numPeriods){
  for(j in i:numPeriods){
    SigmaW[i,j] <- 0.5^(abs(i-j))
    SigmaW[j,i] <- 0.5^(abs(i-j))
  }
}

what <- MASS::mvrnorm(n = 1,
                      mu = rep(mu, numPeriods),
                      Sigma = diag(numPeriods) + SigmaW )

mle <- eb_params_ar1(w=what,SigmaW = SigmaW)

expect_equal(mle$mu,mu, tol = 5/sqrt(numPeriods))
expect_equal(mle$sigmasq,1, tol = 5/sqrt(numPeriods))


# Test that VCV for AR(1) with rho=0.999 close to VCV for AR(1) with rho = 1
expect_equal(
create_V_AR1(rho = 1,
             sigma = 1,
             tVec = c(-10:-1,1:10),
             referencePeriod = 0),
create_V_AR1(rho = 0.999,
             sigma = 1,
             tVec = c(-10:-1,1:10),
             referencePeriod = 0),
tol = 0.1)