library(testthat)
library(MASS)
library(here)
source(here("create_V_AR1.R"))
Tsims <- 100
numSims <- 10^6
rho <- 0.5

eps_mat <-
MASS::mvrnorm(n = numSims,
             Sigma = diag(Tsims),
             mu = rep(0,Tsims))

ar1_mat <-
  matrix(0,
         nrow = NROW(eps_mat),
         ncol = NCOL(eps_mat))

for(t in 1:Tsims){
  rhovec <- matrix(rho^(seq(from = t-1, to =0)),ncol =1 )
  ar1_mat[,t] <- as.matrix(eps_mat[,1:t]) %*% rhovec 
}

basePeriod <- Tsims - 10
delta_mat <- ar1_mat[,c((basePeriod - 10):(basePeriod-1), (basePeriod+1):Tsims )]-
             ar1_mat[,basePeriod]

covDelta <- cov(delta_mat)

V <- create_V_AR1(rho = rho,
                  sigma = 1,
                  tVec = c(-10:-1,1:10))

expect_equal(
  max(abs(V/covDelta-1)),
  0,
  0.1 #allow for 10% sim error in covariance
)
