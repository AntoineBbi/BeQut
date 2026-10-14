#' \code{lslqmm} fits location-scale linear quantile mixed model
#'
#' Function using 'JAGS' software to estimate the linear quantile mixed model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without left formula part
#' @param formFixedScale formula for fixed part of scale regression part without left formula part
#' @param formRandomScale formula for random part of scale regression part without left formula part
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param data long format dataset of observed variables
#' @param tau the quantile's order to be estimated. This must be a number between 0 and 1. If tau is NULL, then the parameter tau is estimated.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param object_lqmm Object returns by 'lqmm' function of BeQuT package to initialize parameter chains, default is \code{NULL}
#' @param precision precision (e.g. inverse of variance) of prior distributions, default is 0.1
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{Blslqmm} object is a list with the following elements:
#'  \describe{
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{jags}}{list including BUGS code and data list used by Jags, and MCMC.list returns by \code{jagsUI::jags()}}
#'   \item{\code{summary_table}}{Data frame returns by \code{MCMCvis::MCMCsummary()} function summarizing posterior statistics of parameter chains.}
#'   \item{\code{out_jagsUI}}{Object returns by \code{jagsUI::jags()} function.}
#'  }
#'
#' @author Antoine Barbieri
#'
#' @import jagsUI
#' @importFrom  MCMCvis MCMCsummary
#'
#' @references Barbieri A, Alcaraz A, Abed M, De Courson H, Jacqmin-Gadda H (2025).
#' \emph{ Asymmetric Laplace distribution regression model for fitting heterogeneous longitudinal response}.
#' arXiv:2512.12362. 2025. doi: https://doi.org/10.48550/arXiv.2512.12362.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' #---- Load data
#' data(threeC)
#'
#' #---- data management
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$age.final65 <- (threeC$age.final-65)/10
#' threeC$age0_65 <- (threeC$age0-65)/10
#' threeC$age.last65 <- (threeC$age.last-65)/10
#' threeC$age.first65 <- (threeC$age.first-65)/10
#' threeC$SBP <- threeC$SBP/10
#'
#' #---- Fit regression model for the median
#' fit_lslqmm <- lslqmm(formFixed = SBP ~ age.visit65,
#'                         formRandom = ~ age.visit65,
#'                         formFixedScale = ~ age.visit65,
#'                         formRandomScale = ~ age.visit65,
#'                         formGroup = ~ ID,
#'                         data = threeC,
#'                         tau = 0.5,
#'                         n.iter = 6000,
#'                         n.burnin = 1000)
#'
#' #---- Get the posterior statistics
#' fit_lslqmm$summary_table
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(fit_lslqmm$jags$out_jags, parameters = "beta")
#'
#' #---- Summary of output
#' summary(fit_lslqmm)
#' }
#'
lslqmm <- function(formFixed,
                   formRandom,
                   formFixedScale,
                   formRandomScale,
                   formGroup,
                   data,
                   tau,
                   n.chains = 3,
                   n.iter = 10000,
                   n.burnin = 5000,
                   n.thin = 1,
                   n.adapt = NULL,
                   object_lqmm = NULL,
                   precision = 0.1,
                   parallel = FALSE){


  # checks of arguments
  if(!is.null(object_lqmm) && !inherits(object_lqmm, "Blqmm"))
    stop("The 'object_lqmm' obejct must be an \"Blqmm\" objects.")
  if(!is.null(object_lqmm) && formFixed!=object_lqmm$control$formFixed)
    stop("The 'formFixed' formula does not match the one from 'object_lqmm'.")

  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),
                             all.vars(formFixed),
                             all.vars(formRandom),
                             all.vars(formFixedScale),
                             all.vars(formRandomScale)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- stats::model.frame(formFixed, data = data_long)
  X1 <- stats::model.matrix(formFixed, mfX)
  mfZ <- stats::model.frame(formRandom, data = data_long)
  Z1 <- stats::model.matrix(formRandom, mfZ)
  mfX <- stats::model.frame(formFixedScale, data = data_long)
  X2 <- stats::model.matrix(formFixedScale, mfX)
  mfZ <- stats::model.frame(formRandomScale, data = data_long)
  Z2 <- stats::model.matrix(formRandomScale, mfZ)
  ncX1 <- ncol(X1)
  ncZ1 <- ncol(Z1)
  ncX2 <- ncol(X2)
  ncZ2 <- ncol(Z2)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))

  # Initialization of prior parameter
  if(!is.null(object_lqmm)){
    # prior beta parameters
    priorMean.beta <- as.vector(object_lqmm$mean$beta)
    priorTau.beta <- if(ncX1==1) c(precision) else diag((1/object_lqmm$StDev$beta^2)*precision)
    # prior xi parameters
    priorMean.xi <- if(ncX2==1) c(log(object_lqmm$mean$sigma)) else c(log(object_lqmm$mean$sigma), rep(0, ncX2-1))
    priorTau.xi <- if(ncX2==1) c(1/sd(log(object_lqmm$sims.list$sigma))*precision) else diag(c(1/sd(log(object_lqmm$sims.list$sigma))*precision,
                                                                                               rep(precision, (length(priorMean.xi)-1) )))
    # priorTau.xi <- if(ncX2==1) c(precision) else diag(rep(precision, length(priorMean.xi)))
  }else{
    # prior beta parameters
    priorMean.beta <- rep(0, ncX1)
    priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
    # prior xi parameters
    priorMean.xi <- if(ncX2==1) c(0) else rep(0, ncX2)
    priorTau.xi <- if(ncX2==1) c(precision) else diag(rep(precision, length(priorMean.xi)))
  }
  # initialization of chains
  initial.values <- list(b = if(is.null(object_lqmm)) matrix(0, ncol = ncZ1, nrow = I) else object_lqmm$random_effect$postMeans,
                         u = matrix(0, ncol = ncZ2, nrow = I),
                         beta = priorMean.beta,
                         xi = priorMean.xi
  )

  # list of data jags
  jags.data <- list(y = y,
                    X1 = X1,
                    Z1 = Z1,
                    X2 = X2,
                    Z2 = Z2,
                    tau = tau,
                    ncX1 = ncX1,
                    ncX2 = ncX2,
                    ncZ1 = ncZ1,
                    ncZ2 = ncZ2,
                    I = I,
                    offset = offset,
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = as.matrix(priorTau.beta),
                    priorMean.xi = priorMean.xi,
                    priorTau.xi = as.matrix(priorTau.xi)
                    )

  if(ncZ1==1 && ncZ2==1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.b = precision,
                        priorB.b = precision,
                        mu_b = 0,
                        priorA.u = precision,
                        priorB.u = precision,
                        mu_u = 0
                   )
    )
    # update initialisation values
    initial.values$prec.Sigma_b <- ifelse(is.null(object_lqmm),
                                          precision,
                                          1/object_lqmm$mean$covariance.b)
    initial.values$prec.Sigma_u <- precision
    # define the appropriate BUGS model
    jags_code <- "model{
  # constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2*va1[j])
    }#end of j loop
    # random effects
    b[i, 1] ~ dnorm(mu_b, prec.Sigma_b)
    u[i, 1] ~ dnorm(mu_u, prec.Sigma_u)
  }#end of i loop
  # priors for parameters
  prec.Sigma_b ~ dgamma(priorA.b, priorB.b)
  covariance.b <- 1/prec.Sigma_b
  prec.Sigma_u ~ dgamma(priorA.u, priorB.u)
  covariance.u <- 1/prec.Sigma_u
  for(p1 in 1:ncX1){
    beta[p1] ~ dnorm(priorMean.beta[p1], priorTau.beta[p1, p1])
  }
  for(p2 in 1:ncX2){
    xi[p2] ~ dnorm(priorMean.xi[p2], priorTau.xi[p2, p2])
  }
}"
  }
  if(ncZ1==1 && ncZ2>1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.b = precision,
                        priorB.b = precision,
                        mu_b = 0,
                        priorA.u = diag(rep(precision, ncZ2)),
                        priorB.u = ncZ2,
                        mu_u = rep(0, ncZ2)
                   )
    )
    # update initialisation values
    if(is.null(object_lqmm)){
      initial.values$prec.Sigma_b <- precision
    }else{
      initial.values$prec.Sigma_b <- 1/object_lqmm$mean$covariance.b
    }
    initial.values$prec.Sigma_u <- diag(precision, ncZ2)
    # define the appropriate BUGS model
    jags_code <- "model{
  # constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2*va1[j])
    }#end of j loop
    # random effects
    b[i, 1] ~ dnorm(mu_b, prec.Sigma_b)
    u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
  }#end of i loop
  # priors for parameters
  prec.Sigma_b ~ dgamma(priorA.b, priorB.b)
  covariance.b <- 1/prec.Sigma_b
  prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)
  covariance.u <- inverse(prec.Sigma_u[, ])
  beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
}"
  }
  if(ncZ1>1 && ncZ2==1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.u = precision,
                        priorB.u = precision,
                        mu_u = 0,
                        priorA.b = diag(rep(precision, ncZ1)),
                        priorB.b = ncZ1,
                        mu_b = rep(0, ncZ1)
                   )
    )
    # update initialisation values
    if(is.null(object_lqmm)){
      initial.values$prec.Sigma_b <- diag(precision, ncol(Z1))
    }else{
      initial.values$prec.Sigma_b <- solve(object_lqmm$mean$covariance.b)
    }
    initial.values$prec.Sigma_u <- precision
    # define the appropriate BUGS model
    jags_code <- "model{
  # constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2*va1[j])
    }#end of j loop
    # random effects
    b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])
    u[i, 1] ~ dnorm(mu_u, prec.Sigma_u)
  }#end of i loop
  # priors for parameters
  prec.Sigma_b[1:ncZ1, 1:ncZ1] ~ dwish(priorA.b[, ], priorB.b)
  covariance.b <- inverse(prec.Sigma_b[, ])
  prec.Sigma_u ~ dgamma(priorA.u, priorB.u)
  covariance.u <- 1/prec.Sigma_u
  beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
}"
  }
  if(ncZ1>1 && ncZ2>1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.b = diag(rep(precision, ncZ1)),
                        priorB.b = ncZ1,
                        mu_b = rep(0, ncZ1),
                        priorA.u = diag(rep(precision, ncZ2)),
                        priorB.u = ncZ2,
                        mu_u = rep(0, ncZ2)
                   )
    )
    # update initialisation values
    if(is.null(object_lqmm)){
      initial.values$prec.Sigma_b <- diag(precision, ncol(Z1))
    }else{
      initial.values$prec.Sigma_b <- solve(object_lqmm$mean$covariance.b)
    }
    initial.values$prec.Sigma_u <- diag(precision, ncol(Z2))
    # define the appropriate BUGS model
    jags_code <- "model{
  # constants
  c1 <- (1-2*tau)/(tau*(1-tau))
  c2 <- 2/(tau*(1-tau))
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2*va1[j])
    }#end of j loop
    # random effects
    b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])
    u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
  }#end of i loop
  # priors for parameters
  prec.Sigma_b[1:ncZ1, 1:ncZ1] ~ dwish(priorA.b[, ], priorB.b)
  covariance.b <- inverse(prec.Sigma_b[, ])
  prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)
  covariance.u <- inverse(prec.Sigma_u[, ])
  beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
}"
  }


  # regression from X1 with beta
  rplc <- paste(paste("beta[", 1:jags.data$ncX1, "] * X1[j, ", 1:jags.data$ncX1, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX1], X1[j, 1:ncX1])", rplc, jags_code, fixed = TRUE)
  # regression on random effects b
  rplc <- paste(paste("b[i, ", 1:jags.data$ncZ1, "] * Z1[j, ", 1:jags.data$ncZ1, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])", rplc, jags_code, fixed = TRUE)
  # regression from X2 with xi
  rplc <- paste(paste("xi[", 1:jags.data$ncX2, "] * X2[j, ", 1:jags.data$ncX2, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(xi[1:ncX2], X2[j, 1:ncX2])", rplc, jags_code, fixed = TRUE)
  # regression on random effects b
  rplc <- paste(paste("u[i, ", 1:jags.data$ncZ2, "] * Z2[j, ", 1:jags.data$ncZ2, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])", rplc, jags_code, fixed = TRUE)

  if(ncZ1==1)
    jags.data$ncZ1 <- NULL
  if(ncZ2==1)
    jags.data$ncZ2 <- NULL

  # initialisation values
  if(n.chains==3)
    inits <- list(initial.values,
                  initial.values,
                  initial.values)
  if(n.chains==2)
    inits <- list(initial.values,
                  initial.values)
  if(n.chains==1)
    inits <- initial.values


  # parameters to save in the sampling step
  parms_to_save <- c("beta", "xi", "b", "covariance.b","u", "covariance.u")

  #---- use JAGS sampler via jagsUI
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = textConnection(jags_code),
                          inits = inits,
                          n.chains = n.chains,
                          parallel = parallel,
                          n.adapt = n.adapt,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          DIC = T)

  #---- output
  out <- list(data = data)
  out$control <- list(formFixed = formFixed,
                      formRandom = formRandom,
                      formGroup = formGroup,
                      formFixedScale = formFixedScale,
                      formRandomScale = formRandomScale,
                      tau = tau,
                      call_function = "lslqmm",
                      n.chains = n.chains,
                      parallel = parallel,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      I = I)

  # save jags information
  out$jags <- list(jags_code = jags_code,
                   jags.data = jags.data,
                   out_jags = out_jags)


  # summary output
  parms_to_save2 <- c("beta", "xi",
                      "covariance.b", "covariance.u")
  out$summary_table <- MCMCvis::MCMCsummary(out$jags$out_jags$samples,
                                            params = parms_to_save2)


  #---- End of the function defining the class and retruning the output
  class(out) <-  c("Blslqmm","Bald", "Bdr", "Bmm")
  out

}
