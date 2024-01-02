library(dplyr)
library(here)
library(ggplot2)

source(here('taupost_posterior.R'))
source(here('create_V_AR1.R'))
source(here('eb_params_ar1.R'))

#Load Benzarti Carloni results
lw_women <- readRDS(here("Data/EmploymentFemale_RefPeriodMinusTwo.rds"))

#Multiply beta and Sigma to convert to pp (i.e. beta * 100,  sigma * 100^2)
beta <- lw_women$beta * 100
sigma <- lw_women$sigma * 100^2
tVec <- lw_women$tVec
referencePeriod <- lw_women$referencePeriod
prePeriodIndices <- which(tVec < referencePeriod)

w <- calc_w(beta = beta,
            Sigma = sigma,
            prePeriodIndices = prePeriodIndices)$w

SigmaW <- calc_w(beta = beta,
            Sigma = sigma,
            prePeriodIndices = prePeriodIndices)$SigmaW
eb_params <- eb_params_ar1(w = w,
                           SigmaW = SigmaW)
