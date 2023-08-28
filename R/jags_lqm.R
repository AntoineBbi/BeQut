jags_lqm <-
  function (...) {
  # constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # quantile linear regression
  for (i in 1:I){
    y[i] ~ dnorm(mu[i], prec[i])
    va1[i] ~ dexp(1/sigma)
    prec[i] <- 1/(sigma*c2*va1[i])
    mu[i] <- inprod(beta[1:ncX], X[i, 1:ncX]) + c1*va1[i]
  }#end of i loop
  # priors for parameters
  for(p in 1:ncX){
    beta[p] ~ dnorm(0, 0.001)
  }
  sigma ~ dgamma(0.001, 0.001)
  }
