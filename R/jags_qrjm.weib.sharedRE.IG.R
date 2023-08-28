jags_qrjm.weib.sharedRE.IG <-
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
    for(r in 1:ncU){
      b[i,r] ~ dnorm(0, prec.Sigma2[r])
    }
    # survival part
    etaBaseline[i] <- inprod(alpha[1:ncZ], Z[i, 1:ncZ]) + inprod(alpha.assoc[1:ncU], b[i, 1:ncU])
    log_S1[i] <- -exp(etaBaseline[i]) * pow(Time[i], shape)
    log_h1[i] <- log(shape) + (shape-1)*log(Time[i]) + etaBaseline[i]
    logL[i] <- event[i]*log_h1[i] + log_S1[i]
    mlogL[i] <- -logL[i] + C
    zeros[i] ~ dpois(mlogL[i])
  }#end of i loop
  # priors for longitudinal parameters
  for(rr in 1:ncU){
    prec.Sigma2[rr] ~ dgamma(priorA.Sigma2, priorB.Sigma2)
    covariance.b[rr] <- 1/prec.Sigma2[rr]
  }
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
  # priors for survival parameters
  alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
  shape ~ dgamma(priorA.shape, priorB.shape)
  for(rr in 1:ncU){
    alpha.assoc[rr] ~ dnorm(0, priorTau.alphaA)
  }
}
