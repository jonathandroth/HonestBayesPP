library(dplyr)
library(here)
library(ggplot2)

source(here('taupost_posterior.R'))
source(here('create_V_AR1.R'))
source(here('eb_params_ar1.R'))

#Load Benzarti Carloni results

lw_women <- readRDS(here("Data/EmploymentFemale_RefPeriodMinusTwo.rds"))
lw_men <- readRDS(here("Data/EmploymentMale_RefPeriodMinusTwo.rds"))

lw_results_fn <-
function(lw_object){
#Multiply beta and Sigma to convert to pp (i.e. beta * 100,  sigma * 100^2)
beta <- lw_object$beta * 100
sigma <- lw_object$sigma * 100^2
tVec <- lw_object$tVec
referencePeriod <- lw_object$referencePeriod
prePeriodIndices <- which(tVec < referencePeriod)

w <- calc_w(beta = beta,
            Sigma = sigma,
            prePeriodIndices = prePeriodIndices)$w

SigmaW <- calc_w(beta = beta,
            Sigma = sigma,
            prePeriodIndices = prePeriodIndices)$SigmaW
eb_params <- eb_params_ar1(w = w,
                           SigmaW = SigmaW)

V_delta_EB <- 
create_V_AR1(rho = 1,
             sigma = sqrt(eb_params$sigmasq),
             tVec = tVec,
             referencePeriod = referencePeriod)

mu_delta_EB <- eb_params$mu * (tVec - referencePeriod)

result <-
taupost_posterior(beta = beta,
                  Sigma = sigma,
                  Vdelta = V_delta_EB,
                  mudelta = mu_delta_EB,
                  tVec = tVec, 
                  Vtaupost = NULL,
                  referencePeriod = referencePeriod)

return(result)

}

lw_results_fn(lw_men)$eventPlot +
  xlim(-10,10) + ylim(-7,3)

lw_results_fn(lw_women)$eventPlot +
  xlim(-10,10) + ylim(-3,6)