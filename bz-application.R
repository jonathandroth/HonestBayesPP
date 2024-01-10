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

posteriorResults <-
taupost_posterior(beta = beta,
                  Sigma = sigma,
                  Vdelta = V,
                  tVec = tVec,
                  referencePeriod = referencePeriod)


eventPlot_withPosterior <-
  posteriorResults$eventPlot +
  ylim(-0.2,0.4) +
  scale_x_continuous(breaks = -4:4) +
  fte_theme()+
  ggplot2::theme(legend.position = c(0.25,0.75))

ggsave(here("Figures/bz-with-posterior.png"),
       width = 6, height = 4)

eventPlot_withoutPosterior <-  
  eventPlot_withPosterior + 
  ggplot2::scale_color_manual(values = c(Original = "#D95F02", Posterior = "transparent"), #make posterior transparent
                              name = NULL)  # Remove the legend title 

  
ggsave(here("Figures/bz-without-posterior.png"),
       width = 6, height = 4)


## The code below calculates the posterior in two ways:
# posterior1 uses an 'informative' prior for taupost with very large variance
# posterior2 uses an 'uninformatve' prior

posterior1 <- taupost_posterior(beta = beta,
                                Sigma = sigma,
                                Vdelta = V,
                                Vtaupost = 10^12 * diag(4),
                                tVec = tVec,
                                referencePeriod = referencePeriod)


posterior2 <- taupost_posterior(beta = beta,
                                Sigma = sigma,
                                Vdelta = V,
                                tVec = tVec, Vtaupost = NULL,
                                referencePeriod = referencePeriod)

posterior1$summaryTable
posterior2$summaryTable


