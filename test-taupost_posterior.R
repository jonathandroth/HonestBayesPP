library(testthat)

# Posterior looks like CIs if deltapre prior is uninformative, deltapost prior variance is zero
result <-
taupost_posterior(beta = c(1,2,3,4),
                  Sigma = diag(4) + 0.5,
                  Vdelta = diag(c(10^8,10^8,10^-8,10^-8) + 10^-10),
                  tVec = c(-2,-1,1,2))

expect_equal(result$summaryTable[,"originalEstimate"],
             result$summaryTable[,"posteriorMean"],
             tol = 10^-6)

expect_equal(result$summaryTable[,"originalSE"],
             result$summaryTable[,"posteriorSD"],
             10^-6)

#As sigma is made very small:
#Posterior variance looks like conditional variance
#Posterior mean looks like betapost - gamma_V * betapre

result <-
taupost_posterior(beta = c(1,2),
                  Sigma = 10^(-8)*diag(2),
                  Vdelta = matrix(c(1,0.5,0.5,1), nrow = 2, ncol =2),
                  tVec = c(-1,1))

expect_equal(sqrt(0.75),
             result$summaryTable[1,"posteriorSD"],
             tol = 10^-6)

expect_equal((2 - 1*0.5),
             result$summaryTable[1,"posteriorMean"],
             tol = 10^-6)



## Run a basic Monte Carlo calibrated to BZ and check coverage of posterior CS ## 

library(dplyr)
library(here)
library(ggplot2)

source(here('taupost_posterior.R'))
source(here('create_V_AR1.R'))
source(here('fte_theme.R'))

#Load Benzarti Carloni results
bz <- readRDS(here("Data/resultsObject-BenzartiCarloni-profits.rds"))

beta <- bz$beta
sigma <- bz$sigma
tVec <- bz$tVec
referencePeriod <- 2008

rho <- 0.766
V <- create_V_AR1(rho = rho,
                  sigma = sqrt(1-rho^2)* 0.063,
                  tVec = tVec)


numSims <- 10^4
delta <- MASS::mvrnorm(n = numSims,Sigma = V, mu = rep(0,nrow(V)))

betahat_error <- MASS::mvrnorm(n = numSims,
                               Sigma = sigma,
                               mu = rep(0,nrow(sigma)))

betahat <- delta + betahat_error
taupost <- 5
betahat[,which(tVec > referencePeriod)] <- betahat[,which(tVec > referencePeriod)] + taupost

simResultsList <-
  purrr::map(.x = 1:numSims,
             .f = ~taupost_posterior(beta = betahat[.x,],
                                     Sigma = sigma,
                                     Vdelta = V,
                                     tVec = tVec,
                                     referencePeriod = referencePeriod))


coverage <- 
  mean(
    purrr::map_dbl(.x = 1:numSims,
                   .f = ~abs((simResultsList[[.x]]$summaryTable$posteriorMean[2]-taupost)/
                               simResultsList[[.x]]$summaryTable$posteriorSD[2])<1.96 )
  )

expect_equal(coverage,
             0.95,
             tol = 0.01)






