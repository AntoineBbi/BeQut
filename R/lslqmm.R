#' \code{lslqmm} fits location-scale linear quantile mixed model
#'
#' Function using 'JAGS' software to estimate the linear quantile mixed model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without left formula part
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param formFixedScale formula for fixed part of scale regression part without left formula part
#' @param formRandomScale formula for random part of scale regression part without left formula part
#' @param data long format dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param precision variance by default for vague prior distribution
#' @param save_jagsUI If \code{TRUE} (by default), the output of \code{jagsUI} package is returned by the function. Warning, if \code{TRUE}, the output can be large.
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{Blqmm} object is a list with the following elements:
#'  \describe{
#'   \item{\code{mean}}{list of posterior mean for each parameter}
#'   \item{\code{median}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95 for each parameters excepted for covariance parameters in covariance matrix of random effects. Otherwise, use save_jagsUI=TRUE to have the associated quantiles.}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters and random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{random_effect}}{list for each quantile including both posterior mean and posterior standard deviation of subject-specific random effects}
#'   \item{\code{out_jagsUI}}{only if \code{save_jagsUI=TRUE} in argument: list including posterior mean, median, quantiles (2.5%, 25%, 50%, 75%, 97.5%), standart deviation for each parameter and each random effect.
#'   Moreover, this list also returns the MCMC draws, the Gelman and Rubin diagnostics (see output of jagsUI objects)}
#'  }
#'
#' @author Antoine Barbieri
#'
#' @import lqmm jagsUI
#'
#' @references Courcoul L, Tzourio C, Woodward M, Barbieri A, Jacqmin-Gadda H (2024).
#' \emph{	A location-scale joint model for studying the link between the time-dependent subject-specific variability of blood pressure and competing events}.
#' arXiv, arXiv:2306.16785. doi: https://doi.org/10.48550/arXiv.2306.16785.
#' @references Marco Geraci and Matteo Bottai (2014).
#' \emph{Linear quantile mixed models}.
#' Statistics and Computing, 24(3):461-479. doi: 10.1007/s11222-013-9381-9.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' #---- Use dataLong dataset
#' data(pbc2, package = "JMbayes")
#'
#' #---- Fit regression model for the median
#' lslqmm_median <- lslqmm(formFixed = log(serBilir) ~ year,
#'                         formRandom = ~ year,
#'                         formFixedScale = ~ year,
#'                         formRandomScale = ~ year,
#'                         formGroup = ~ id,
#'                         data = pbc2,
#'                         tau = 0.5,
#'                         n.iter = 10000,
#'                         n.burnin = 1000)
#'
#' #---- Get the posterior means
#' lslqmm_median$mean
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(lslqmm_median$out_jagsUI, parameters = "beta")
#'
#' #---- Summary of output
#' summary(lslqmm_median)
#' }
#'
lslqmm <- function(formFixed,
                   formRandom,
                   formGroup,
                   formFixedScale,
                   formRandomScale,
                   data,
                   tau,
                   n.chains = 3,
                   n.iter = 10000,
                   n.burnin = 5000,
                   n.thin = 1,
                   n.adapt = NULL,
                   precision = 10,
                   save_jagsUI = TRUE,
                   parallel = FALSE){


  #-- To do list
  # change the example
  # iniitalisation parameters using BeQut

  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- stats::model.frame(formFixed, data = data_long)
  X1 <- stats::model.matrix(formFixed, mfX)
  mfZ <- stats::model.frame(formRandom, data = data_long)
  Z1 <- stats::model.matrix(formRandom, mfZ)
  mfX <- stats::model.frame(formFixedScale, data = data_long)
  X2 <- stats::model.matrix(formFixedScale, mfX)
  mfZ <- stats::model.frame(formRandomScale, data = data_long)
  Z2 <- stats::model.matrix(formRandomScale, mfZ)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  # use lqmm function to initiated values
  message("Initiation of parameter values using lqmm package.")
  tmp_model <- lqmm::lqmm(fixed = formFixed,
                          random = formRandom,
                          group = id,
                          tau = tau,
                          data = data_long)
  # prior beta parameters
  priorMean.beta <- lqmm::coef.lqmm(tmp_model)
  priorTau.beta <- diag(rep(1/precision, length(priorMean.beta)))
  # prior xi parameters
  ncX2 <- ncol(X2)
  priorMean.xi <- if(ncX2==1) tmp_model$scale else c(log(tmp_model$scale), rep(0, ncX2-1))
  priorTau.xi <- diag(rep(1/10, length(priorMean.beta)))

  bis <- as.matrix(lqmm::ranef(tmp_model))
  bis[abs(bis)<.0001] <- 0
  uis <- matrix(0, ncol = ncX2, nrow = I)
  initial.values <- list(b = bis,
                         u = uis,
                         beta = priorMean.beta,
                         xi = priorMean.xi)

  # list of data jags
  jags.data <- list(y = y,
                    X1 = X1,
                    Z1 = Z1,
                    X2 = X2,
                    Z2 = Z2,
                    tau = tau,
                    ncX1 = ncol(X1),
                    ncZ1 = ncol(Z1),
                    ncX2 = ncol(X2),
                    ncZ2 = ncol(Z2),
                    I = I,
                    offset = offset,
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = priorTau.beta,
                    priorMean.xi = priorMean.xi,
                    priorTau.xi = priorTau.xi
  )

  if(jags.data$ncZ1==1 && jags.data$ncZ2==1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.b = 1/precision,
                        priorB.b = 1/precision,
                        mu_b = 0,
                        priorA.u = 1/precision,
                        priorB.u = 1/precision,
                        mu_u = 0
                   )
    )
    # update initialisation values
    initial.values$prec.Sigma_b <- 1/lqmm::VarCorr(tmp_model)
    initial.values$prec.Sigma_u <- 1/precision
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
  beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
}"
  }
  if(jags.data$ncZ1==1 && jags.data$ncZ2>1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.b = 1/precision,
                        priorB.b = 1/precision,
                        mu_b = 0,
                        priorA.u = diag(rep(1/precision, ncol(Z2))),
                        priorB.u = ncol(Z2),
                        mu_u = rep(0, ncol(Z2))
                   )
    )
    # update initialisation values
    initial.values$prec.Sigma_b <- 1/lqmm::VarCorr(tmp_model)
    initial.values$prec.Sigma_u <- diag(1/precision, ncol(Z2))
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
  if(jags.data$ncZ1>1 && jags.data$ncZ2==1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.u = 1/precision,
                        priorB.u = 1/precision,
                        mu_u = 0,
                        priorA.b = diag(rep(1/precision, ncol(Z1))),
                        priorB.b = ncol(Z1),
                        mu_b = rep(0, ncol(Z1))
                        )
                   )
    # update initialisation values
    initial.values$prec.Sigma_b <- diag(1/lqmm::VarCorr(tmp_model))
    initial.values$prec.Sigma_u <- 1/precision
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
  if(jags.data$ncZ1>1 && jags.data$ncZ2>1){
    # update jags.data list
    jags.data <- c(jags.data,
                   list(priorA.b = diag(rep(1/precision, ncol(Z1))),
                        priorB.b = ncol(Z1),
                        mu_b = rep(0, ncol(Z1)),
                        priorA.u = diag(rep(1/precision, ncol(Z2))),
                        priorB.u = ncol(Z2),
                        mu_u = rep(0, ncol(Z2))
                        )
                   )
    # update initialisation values
    initial.values$prec.Sigma_b <- diag(1/lqmm::VarCorr(tmp_model))
    initial.values$prec.Sigma_u <- diag(1/precision, ncol(Z2))
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
                      tau = tau,
                      call_function = "lqmm",
                      n.chains = n.chains,
                      parallel = parallel,
                      n.adapt = n.adapt,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      I = I)

  #- other outputs

  # sims.list output
  out$sims.list <- out_jags$sims.list
  out$sims.list$b <- NULL
  out$sims.list$u <- NULL

  # b random effect output
  out$random_effect_b <- list(postMeans = out_jags$mean$b,
                              postSd = out_jags$sd$b,
                              Rhat = out_jags$Rhat$b)
  colnames(out$random_effect_b$postMeans) <-
    colnames(out$random_effect_b$postSd) <-
    colnames(out$random_effect_b$Rhat) <- colnames(Z1)

  # u random effect output
  out$random_effect_u <- list(postMeans = out_jags$mean$u,
                              postSd = out_jags$sd$u,
                              Rhat = out_jags$Rhat$u)
  colnames(out$random_effect_u$postMeans) <-
    colnames(out$random_effect_u$postSd) <-
    colnames(out$random_effect_u$Rhat) <- colnames(Z2)

  # median : Posterior median of parameters
  out$median <- out_jags$q50
  out$median$b <- NULL
  out$median$u <- NULL

  # mean : Posterior mean of parameters
  out$mean <- out_jags$mean
  out$mean$b <- NULL
  out$mean$u <- NULL

  # modes of parameters
  out$modes <- lapply(out$sims.list, function(x) {
    m <- function(x) {
      d <- stats::density(x, bw = "nrd", adjust = 3, n = 1000)
      d$x[which.max(d$y)]
    }
    if (is.matrix(x))
      as.array(apply(x, 2, m))
    else{
      if(is.array(x))
        apply(x, c(2,3), m)
      else m(x)
    }
  })

  # standard error of parameters
  out$StErr <- lapply(out$sims.list, function(x) {
    f <- function(x) {
      acf.x <- drop(stats::acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
      acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
      ess <- length(x)/(1 + 2 * sum(acf.x))
      sqrt(stats::var(x)/ess)
    }
    if (is.matrix(x))
      as.array(apply(x, 2, f))
    else{
      if(is.array(x))
        apply(x, c(2,3), f)
      else f(x)
    }
  })

  # standard deviation of parameters
  out$StDev <- out_jags$sd
  out$StDev$b <- NULL
  out$StDev$u <- NULL

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  out$Rhat$b <- NULL
  out$Rhat$u <- NULL

  # names
  names(out$mean$beta) <-
    names(out$median$beta) <-
    names(out$modes$beta) <-
    names(out$StErr$beta) <-
    names(out$Rhat$beta) <-
    names(out$StDev$beta) <- colnames(X1)

  names(out$mean$xi) <-
    names(out$median$xi) <-
    names(out$modes$xi) <-
    names(out$StErr$xi) <-
    names(out$Rhat$xi) <-
    names(out$StDev$xi) <- colnames(X2)

  if(ncol(Z1)==1){
    names(out$mean$covariance.b) <-
      names(out$median$covariance.b) <-
      names(out$modes$covariance.b) <-
      names(out$StErr$covariance.b) <-
      names(out$Rhat$covariance.b) <-
      names(out$StDev$covariance.b) <- colnames(Z1)
  }else{
    colnames(out$mean$covariance.b) <-
      rownames(out$mean$covariance.b) <-
      colnames(out$median$covariance.b) <-
      rownames(out$median$covariance.b) <-
      colnames(out$modes$covariance.b) <-
      rownames(out$modes$covariance.b) <-
      colnames(out$Rhat$covariance.b) <-
      rownames(out$Rhat$covariance.b) <-
      colnames(out$StErr$covariance.b) <-
      rownames(out$StErr$covariance.b) <-
      colnames(out$StDev$covariance.b) <-
      rownames(out$StDev$covariance.b) <- colnames(Z1)
  }

  if(ncol(Z2)==1){
    names(out$mean$covariance.u) <-
      names(out$median$covariance.u) <-
      names(out$modes$covariance.u) <-
      names(out$StErr$covariance.u) <-
      names(out$Rhat$covariance.u) <-
      names(out$StDev$covariance.u) <- colnames(Z2)
  }else{
    colnames(out$mean$covariance.u) <-
      rownames(out$mean$covariance.u) <-
      colnames(out$median$covariance.u) <-
      rownames(out$median$covariance.u) <-
      colnames(out$modes$covariance.u) <-
      rownames(out$modes$covariance.u) <-
      colnames(out$Rhat$covariance.u) <-
      rownames(out$Rhat$covariance.u) <-
      colnames(out$StErr$covariance.u) <-
      rownames(out$StErr$covariance.u) <-
      colnames(out$StDev$covariance.u) <-
      rownames(out$StDev$covariance.u) <- colnames(Z2)
  }

  # credible intervalles

  ## beta parameters
  out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                        as.vector(t(out_jags$q97.5$beta)))
  rownames(out$CIs$beta) <- colnames(X1)
  colnames(out$CIs$beta) <- c("2.5%", "97.5%")

  ## xi parameters
  out$CIs$xi <- cbind(as.vector(t(out_jags$q2.5$xi)),
                        as.vector(t(out_jags$q97.5$xi)))
  rownames(out$CIs$xi) <- colnames(X1)
  colnames(out$CIs$xi) <- c("2.5%", "97.5%")

  # only for diagonal elements of covariance matrix of random effects
  if(ncol(Z1)==1){
    out$CIs$covariance.b <- cbind(as.vector(t(out_jags$q2.5$covariance.b)),
                                  as.vector(t(out_jags$q97.5$covariance.b)))
    rownames(out$CIs$covariance.b) <- colnames(Z1)
    colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  }else{
    out$CIs$covariance.b <- cbind(as.vector(diag(out_jags$q2.5$covariance.b)),
                                  as.vector(diag(out_jags$q97.5$covariance.b)))
    rownames(out$CIs$covariance.b) <- colnames(Z1)
    colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  }
  if(ncol(Z2)==1){
    out$CIs$covariance.u <- cbind(as.vector(t(out_jags$q2.5$covariance.u)),
                                  as.vector(t(out_jags$q97.5$covariance.u)))
    rownames(out$CIs$covariance.u) <- colnames(Z2)
    colnames(out$CIs$covariance.u) <- c("2.5%", "97.5%")
  }else{
    out$CIs$covariance.u <- cbind(as.vector(diag(out_jags$q2.5$covariance.u)),
                                  as.vector(diag(out_jags$q97.5$covariance.u)))
    rownames(out$CIs$covariance.u) <- colnames(Z2)
    colnames(out$CIs$covariance.u) <- c("2.5%", "97.5%")
  }

  # save jags output if requires
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  #---- End of the function defining the class and retruning the output
  class(out) <- "Blslqmm"
  out

}
