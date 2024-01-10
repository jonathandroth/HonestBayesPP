library(dplyr)
library(here)
library(ggplot2)

source(here('taupost_posterior.R'))
source(here('create_V_AR1.R'))
source(here('eb_params_ar1.R'))
source(here('fte_theme.R'))

#Load Benzarti Carloni results

lw_women <- readRDS(here("Data/EmploymentFemale_RefPeriodMinusTwo.rds"))
lw_men <- readRDS(here("Data/EmploymentMale_RefPeriodMinusTwo.rds"))

lw_results_fn <-
function(lw_object,
         useMLE = T){
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
                           SigmaW = SigmaW,
                           useMLE = useMLE)

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

eb_params$sigma <- sqrt(eb_params$sigmasq)
result$eb_params <- eb_params

return(result)

}

eventPlot_withPosterior_men <- 
lw_results_fn(lw_men)$eventPlot +
  xlim(-10,10.5) + ylim(-7,3) +
  fte_theme()+
  ggplot2::theme(legend.position = c(0.25,0.25))


ggsave(here("Figures/lz-men-with-posterior.png"),
       width = 6, height = 4)

  
eventPlot_withoutPosterior_men <-   
  eventPlot_withPosterior_men + 
  ggplot2::scale_color_manual(values = c(Original = "#D95F02", Posterior = "transparent"), #make posterior transparent
                              name = NULL)  # Remove the legend title 

ggsave(here("Figures/lz-men-without-posterior.png"),
       width = 6, height = 4)


#Male EB Params
lw_results_fn(lw_men)$eb_params


eventPlot_withPosterior_women <- 
lw_results_fn(lw_women)$eventPlot +
  xlim(-10,10.5) + ylim(-3,6) +
  fte_theme()+
  ggplot2::theme(legend.position = c(0.25,0.75))


ggsave(here("Figures/lz-women-with-posterior.png"),
       width = 6, height = 4)


eventPlot_withoutPosterior_women <-   
  eventPlot_withPosterior_women + 
  ggplot2::scale_color_manual(values = c(Original = "#D95F02", Posterior = "transparent"), #make posterior transparent
                              name = NULL)  # Remove the legend title 

ggsave(here("Figures/lz-women-without-posterior.png"),
       width = 6, height = 4)

#Female EB Params
lw_results_fn(lw_women)$eb_params
