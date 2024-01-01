#' @description
#' This function computes the posterior for taupost under a normal prior
#' @param beta The vector of event-study coefficients
#' @param Sigma The variance-covariance matrix of the event-study coefficients
#' @param Vdelta The prior variance-covariance matrix on delta (the violation of PT)
#' @param Vtaupost (Optional) The prior variance-covariance matrix on the vector of treatment effects. 
#' Null signifies uninformative prior. Default is NULL.
#' @param mutaupost (Optional) The prior mean for taupost. Default is vector of zeros (which shouldn't matter if using uninformative prior)
#' @param mudelta (Optional) The prior mean for delta. Default is a vector of zeros

taupost_posterior <- function(beta, 
                              Sigma, 
                              Vdelta,
                              mudelta = matrix(0,ncol = 1, nrow = length(tVec)),
                              Vtaupost = NULL,
                              mutaupost = matrix(0, nrow =sum(tVec - referencePeriod > 0), ncol = 1),
                              tVec, 
                              referencePeriod = 0){
  
  relativeTime <- tVec - referencePeriod
  prePeriodIndices <- which(relativeTime < 0)
  postPeriodIndices <- which(relativeTime > 0)
  
  numPrePeriods <- length(prePeriodIndices)
  numPostPeriods <- length(postPeriodIndices)
  
  if(!is.null(Vtaupost)){
  
  #Create a vcv for the prior on tau = (0',taupost')
  Vtau <- matrix(0, 
                 nrow = numPrePeriods + numPostPeriods,
                 ncol = numPrePeriods + numPostPeriods)
  Vtau[postPeriodIndices, postPeriodIndices] <- Vtaupost
  
  #Create the vector for the prior mean on tau, mutau = (0',mutaupost)
  mutau <- rbind(matrix(0, nrow = numPrePeriods, ncol = 1),
                 mutaupost)
    
  #Create the prior on beta = delta + tau (these are assumed indep, so we sum the variances)  
  Vbeta <- Vdelta + Vtau
  mubeta <- mudelta + mutau
  
  
  #Posterior mean and variance for beta
  betaPosterior <- solve( solve(Vbeta) + solve(Sigma) ) %*% ( solve(Sigma) %*% beta + solve(Vbeta) %*% mubeta  )
  VbetaPosterior <- solve( solve(Vbeta) + solve(Sigma) )
  
  #Calculate the posterior covariance between beta and taupost
  # We have Cov(betapre, taupost) = 0 and Cov(betapost, taupost) = Vtaupost
  Vbetatau <- matrix(0,
                     nrow = numPostPeriods,
                     ncol = numPostPeriods + numPrePeriods)
  
  Vbetatau[ , postPeriodIndices] <- Vtaupost
  
  #Calculate the posterior for tau
  tauPostPosterior <- mutaupost + Vbetatau %*% solve(Vbeta) %*% (betaPosterior - mubeta)
  VtauPostPosterior <- (Vtaupost - Vbetatau %*% solve(Vbeta) %*% t(Vbetatau) ) +
    Vbetatau %*% solve(Vbeta) %*% VbetaPosterior %*% t(Vbetatau %*% solve(Vbeta))
  
  }else{
    Sigmapre <- Sigma[prePeriodIndices, prePeriodIndices]
    Sigmapost <- Sigma[postPeriodIndices, postPeriodIndices]
    Sigmaprepost <- Sigma[prePeriodIndices, postPeriodIndices]
    GammaSigma <- solve(Sigmapre) %*% Sigmaprepost
    
    Vpre <- Vdelta[prePeriodIndices, prePeriodIndices]
    Vpost <- Vdelta[postPeriodIndices, postPeriodIndices]
    Vprepost <- Vdelta[prePeriodIndices, postPeriodIndices]
    GammaV <- solve(Vpre) %*% Vprepost  
    
    mupre <- mudelta[prePeriodIndices]
    mupost <- mudelta[postPeriodIndices]
    
    betapre <- beta[prePeriodIndices]
    betapost <- beta[postPeriodIndices]
    
    betapreposterior <- solve( solve(Sigmapre) + solve(Vpre)  ) %*% 
                         (solve(Sigmapre) %*% betapre + solve(Vpre) %*% mupre)
    
    betapostposterior <- betapost - t(GammaSigma) %*% (betapre - betapreposterior)
    
    betaposterior <- matrix(NA, nrow = numPrePeriods + numPostPeriods)
    betaposterior[prePeriodIndices] <- betapreposterior
    betaposterior[postPeriodIndices] <- betapostposterior
    
    Vbetapreposterior <- solve( solve(Sigmapre) + solve(Vpre)  ) 
    Vbetapostposterior <- Sigmapost - t(Sigmaprepost) %*% solve(Sigmapre) %*% Sigmaprepost +
                          t(GammaSigma) %*% Vbetapreposterior %*% GammaSigma
    
    Vbetaprepostposterior <- t(GammaSigma) %*% Vbetapreposterior
    
    Vbetaposterior <- matrix(NA,
                             nrow = numPrePeriods + numPostPeriods,
                             ncol = numPrePeriods + numPostPeriods
                             )
    
    Vbetaposterior[prePeriodIndices, prePeriodIndices] <- Vbetapreposterior
    Vbetaposterior[postPeriodIndices, postPeriodIndices] <- Vbetapostposterior
    Vbetaposterior[prePeriodIndices, postPeriodIndices] <- Vbetaprepostposterior
    Vbetaposterior[postPeriodIndices, prePeriodIndices] <- t(Vbetaprepostposterior)
    
    tauPostPosterior <- betapostposterior - mupost - t(GammaV) %*% (betapreposterior - mupre)
    
    #Weights on betapre and betapost that give E[delta | beta]
    Wpost <- diag(numPostPeriods)
    Wpre <- t(GammaV)
    W <- matrix(NA, ncol = numPrePeriods + numPostPeriods, nrow = numPostPeriods)
    W[, prePeriodIndices] <- Wpre
    W[, postPeriodIndices] <- Wpost
    
    VtauPostPosterior <- Vpost - t(Vprepost) %*% solve(Vpre) %*% Vprepost +
                        W %*% Vbetaposterior %*% t(W)
  }
  
  summaryTable <- data.frame(relativeTime = relativeTime[postPeriodIndices],
                             originalEstimate = beta[postPeriodIndices],
                             originalSE = sqrt(diag(Sigma)[postPeriodIndices]),
                             posteriorMean = tauPostPosterior,
                             posteriorSD = sqrt(diag(VtauPostPosterior)))
  
  #Create an event-study plot
  
  #Put original estimates in a df
  eventstudyDF_original <- data.frame(relativeTime = relativeTime, 
                                      Estimate = beta,
                                      SE = sqrt(diag(Sigma)),
                                      type = "Original")
  
  #Add reference period
  eventstudyDF_original <- rbind(eventstudyDF_original,
                                 data.frame(relativeTime = 0, 
                                            Estimate = 0,
                                            SE = NA,
                                            type = "Original"))
  
  #Put posteriors in a df with same columns but different type
  # Note: SE <-> posterior SD
  eventstudyDF_posterior <- data.frame(relativeTime = relativeTime[postPeriodIndices],
                                       Estimate = tauPostPosterior,
                                       SE = sqrt((diag(VtauPostPosterior))),
                                       type = "Posterior")
  
  eventstudyDF <- rbind(eventstudyDF_original,
                        eventstudyDF_posterior)
  

  #Create event-study
  dodge_width <- 0.5  # Adjust this value as needed

  eventPlot <- ggplot2::ggplot(data = eventstudyDF,
                               ggplot2::aes(x = relativeTime,
                                            y = Estimate,
                                            ymin = Estimate - 1.96*SE,
                                            ymax = Estimate + 1.96*SE,
                                            color = type)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = dodge_width)) +
    ggplot2::geom_errorbar(width = 0.1,  # Adjust this value for horizontal bar width
                           position = ggplot2::position_dodge(width = dodge_width)) +
    ggplot2::scale_color_manual(values = c("Original" = "#D95F02", Posterior = "#666666"),
                                name = NULL)  # Remove the legend title
  

  
  # #Create tauposterior = (NA, taupost)
  # tauPosterior <- matrix(NA, nrow = numPrePeriods + numPostPeriods)
  # tauPosterior[postPeriodIndices] <- tauPostPosterior
  # 
  # tauPosteriorSD <- matrix(NA, nrow = numPrePeriods + numPostPeriods)
  # tauPosteriorSD[postPeriodIndices] <- sqrt(diag(VtauPostPosterior))
  # 
  # eventstudyDF <- data.frame(relativeTime = relativeTime, 
  #                            originalEstimate = beta,
  #                            originalSE = sqrt(diag(Sigma)),
  #                            posteriorMean = tauPosterior,
  #                            posteriorSD = tauPosteriorSD)
  # 
  # #Add row with referencePeriod
  # eventstudyDF <- rbind(eventstudyDF,
  #                       data.frame(relativeTime = 0, 
  #                                  originalEstimate = 0,
  #                                  originalSE = NA,
  #                                  posteriorMean = NA,
  #                                  posteriorSD = NA)
  #                 )
  # 
  # eventstudyDF$originalUB <- eventstudyDF$originalEstimate + 1.96 * eventstudyDF$originalSE
  # eventstudyDF$originalLB <- eventstudyDF$originalEstimate - 1.96 * eventstudyDF$originalSE
  # 
  # eventstudyDF$posteriorUB <- eventstudyDF$posteriorMean + 1.96 * eventstudyDF$posteriorSD
  # eventstudyDF$posteriorLB <- eventstudyDF$posteriorMean - 1.96 * eventstudyDF$posteriorSD
  # 
  # 
  # dodge_width <- 0.5  # Adjust this value as needed
  # 
  # eventPlot <-
  # ggplot2::ggplot(data = eventstudyDF) +
  #   ggplot2::geom_point(ggplot2::aes(x = relativeTime, y = originalEstimate, 
  #                                    color = "Original", group = "Original"), 
  #                       position = ggplot2::position_dodge(width = dodge_width)) +
  #   ggplot2::geom_errorbar(ggplot2::aes(x = relativeTime, 
  #                                       ymin = originalLB, ymax = originalUB,
  #                                       color = "Original", group = "Original"), 
  #                          width = 0.2, 
  #                          position = ggplot2::position_dodge(width = dodge_width)) +
  #   ggplot2::geom_point(ggplot2::aes(x = relativeTime, y = posteriorMean, 
  #                                    color = "Posterior", group = "Posterior"), 
  #                       position = ggplot2::position_dodge(width = dodge_width)) +
  #   ggplot2::geom_errorbar(ggplot2::aes(x = relativeTime, 
  #                                       ymin = posteriorLB, ymax = posteriorUB,
  #                                       color = "Posterior", group = "Posterior"), 
  #                          width = 0.2, 
  #                          position = ggplot2::position_dodge(width = dodge_width)) +
  #   ggplot2::scale_color_manual(values = c("Original" = "blue", "Posterior" = "orange"))
  # 
  
  return(list(summaryTable = summaryTable,
              tauPosterior = tauPostPosterior,
              sigmaPosterior = VtauPostPosterior,
              eventPlot = eventPlot,
              eventPlotDF = eventstudyDF))
}