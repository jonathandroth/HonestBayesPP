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

