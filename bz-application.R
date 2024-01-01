library(dplyr)
library(here)
library(ggplot2)

source(here('taupost_posterior.R'))
source(here('create_V_AR1.R'))

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

taupost_posterior(beta = beta,
                  Sigma = sigma,
                  Vdelta = V,
                  tVec = tVec,
                  referencePeriod = referencePeriod)$eventPlot +
  ylim(-0.2,0.4)
  



taupost_posterior(beta = beta,
                  Sigma = sigma,
                  Vdelta = V,
                  tVec = tVec, Vtaupost = NULL,
                  referencePeriod = referencePeriod)$eventPlot +
  ylim(-0.2,0.4)



# posterior1 <- taupost_posterior(beta = beta,
#                                 Sigma = sigma,
#                                 Vdelta = V,
#                                 Vtaupost = 10^12 * diag(4),
#                                 tVec = tVec,
#                                 referencePeriod = referencePeriod)
# 
# 
# posterior2 <- taupost_posterior(beta = beta,
#                                 Sigma = sigma,
#                                 Vdelta = V,
#                                 tVec = tVec, Vtaupost = NULL,
#                                 referencePeriod = referencePeriod)
# 
# posterior1$summaryTable
# posterior2$summaryTable
# 
# 
