jags_qrjm.weib.value <-
  function (...) {
  # constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      va1[j] ~ dexp(1/sigma)
      prec[j] <- 1/(sigma*c2*va1[j])
      mu[j] <- inprod(beta[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncU], U[j, 1:ncU]) + c1*va1[j]
    }#end of j loop
    # random effects
    b[i, 1:ncU] ~ dmnorm(mu0[], prec.Sigma2[, ])
    # survival part
    etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ])
    shareY[i] <- inprod(beta[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 1:ncU], Utime[i, 1:ncU])
    log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i] + alpha.assoc * shareY[i]
    for (k in 1:K) {
      shareY.s[i, k] <- inprod(beta[1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(b[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
      SurvLong[i, k] <- wk[k] * shape * pow(st[i, k], shape - 1) * exp(alpha.assoc * shareY.s[i, k])
    }
    log_S1[i] <- (-exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ]))
    logL[i] <- event[i]*log_h1[i] + log_S1[i]
    mlogL[i] <- -logL[i] + C
    zeros[i] ~ dpois(mlogL[i])
  }#end of i loop
  # priors for longitudinal parameters
  prec.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b <- inverse(prec.Sigma2[, ])
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
  # priors for survival parameters
  alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
  shape ~ dgamma(priorA.shape, priorB.shape)
  alpha.assoc ~ dnorm(0, priorTau.alphaA)
}
