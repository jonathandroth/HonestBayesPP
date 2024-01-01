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
  

