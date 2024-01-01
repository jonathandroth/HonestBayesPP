a <- 10^7
sigmapost <- 1+a
sigmapre <- 3
sigmaprepost <- 0.5

vpre <- 1
vpost <- 2
vprepost <- 0.6


Sigma <- matrix(c(sigmapre,sigmaprepost,sigmaprepost,sigmapost),
                nrow = 2,
                ncol = 2,
                byrow = T)

V <- matrix(c(vpre,vprepost,vprepost,vpost),
            nrow = 2,
            ncol = 2,
            byrow = T)

detSigma <- determinant.matrix(x = Sigma, logarithm = F)$modulus[1]
detV <- determinant.matrix(x = V, logarithm = F)$modulus[1]
detSigmaPlusVInv <- determinant.matrix(x = solve(V) + solve(Sigma), logarithm = F)$modulus[1]

solve(solve(V) + solve(Sigma)) %*% solve(Sigma)
solve(solve(V) + solve(Sigma)) %*% solve(V)

(1/vpre)/(1/sigmapre + 1/vpre) * matrix(c(1,0,-vprepost/sigmapre,1+vpre/sigmapre),
                                        nrow = 2,
                                        byrow = T)


# #Check that we get (Sigma^-1 + V^-1)^-1 right
# solve(solve(V) + solve(Sigma))
# 
# 1/(detSigmaPlusVInv) * matrix(c(1/detSigma * sigmapre + 1/detV * vpre,
#                                 1/detSigma * sigmaprepost + 1/detV * vprepost,
#                                 1/detSigma * sigmaprepost + 1/detV * vprepost,
#                                 1/detSigma * sigmapost + 1/detV * vpost
#                                 ),
#                               nrow = 2, byrow = T)
# 
# #Looks good
# 
# 
# #Check that we get the bottom left coef right without simplifying
# (solve(solve(V) + solve(Sigma)) %*% solve(V))[2,1]
# (1/detV)/detSigmaPlusVInv* ((1/detSigma * sigmaprepost + 1/detV * vprepost)*vpost - (1/detSigma * sigmapost + 1/detV * vpost) * vprepost)
