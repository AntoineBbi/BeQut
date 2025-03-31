#' \code{scALDreg} fits density regression for longitudinal data assuming asymmetric Laplace distribution
#'
#' Function using 'JAGS' software to estimate the individual response distribution assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formFixedScale formula for fixed part of scale regression with left formula part is response variable (log link)
#' @param formRandomScale formula for random part of scale regression without left formula part  (log link)
#' @param Location numeric value, must be'0' or '1' to fix Location part to 0 or to fit an intercept, respectively."
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param tau skewness (vector) parameter. By default, this is a scalar between 0 and 1, but it can be a n-length vector specifying the individual skewness.
#' @param lambda numeric regarding the penanlisation on fixed scale part. By default, 0.
#' @param data long format dataset of observed variables
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param precision precision (e.g. inverse of variance) of prior distributions, default is 0.1
#' @param save_jagsUI If \code{TRUE} (by default), the output of \code{jagsUI} package is returned by the function. Warning, if \code{TRUE}, the output can be large.
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{BpscALDreg} object is a list with the following elements:
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
#'   \item{\code{random_effect_b}}{list for each quantile including both posterior mean and posterior standard deviation of subject-specific location random effects}
#'   \item{\code{random_effect_u}}{list for each quantile including both posterior mean and posterior standard deviation of subject-specific scale random effects}
#'   \item{\code{random_effect_tau}}{list for each quantile including both posterior mean and posterior standard deviation of subject-specific skewness random effects}
#'   \item{\code{out_jagsUI}}{only if \code{save_jagsUI=TRUE} in argument: list including posterior mean, median, quantiles (2.5%, 25%, 50%, 75%, 97.5%), standart deviation for each parameter and each random effect.
#'   Moreover, this list also returns the MCMC draws, the Gelman and Rubin diagnostics (see output of jagsUI objects)}
#'  }
#'
#' @author Antoine Barbieri
#'
#' @import jagsUI
#'
#' @references Courcoul L, Tzourio C, Woodward M, Barbieri A, Jacqmin-Gadda H (2024).
#' \emph{	A location-scale joint model for studying the link between the time-dependent subject-specific variability of blood pressure and competing events}.
#' arXiv, arXiv:2306.16785. doi: https://doi.org/10.48550/arXiv.2306.16785.
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
#' test <- scALDreg(formFixedScale = log(serBilir) ~ year,
#'                   formRandomScale = ~ year,
#'                   formGroup = ~ id,
#'                   data = pbc2,
#'                   tau = 0.5,
#'                   lambda = 0,
#'                   n.iter = 2000,
#'                   n.burnin = 1000)
#'
#' #---- Get the posterior means
#' test$mean
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(test$out_jagsUI, parameters = "xi")
#'
#' #---- Summary of output
#' summary(test)
#' }
#'
scALDreg <- function(formFixedScale,
                    formRandomScale,
                    Location = 0,
                    formGroup,
                    tau,
                    lambda = 0,
                    data,
                    n.chains = 3,
                    n.iter = 6000,
                    n.burnin = 1000,
                    n.thin = 1,
                    n.adapt = NULL,
                    precision = 0.1,
                    save_jagsUI = TRUE,
                    parallel = FALSE){
  
  
  # #---- Load data
  # data(threeC, package = "LSJM")
  # 
  # #---- data management
  # threeC$age.visit65 <- (threeC$age.visit-65)/10
  # threeC$age.final65 <- (threeC$age.final-65)/10
  # threeC$age0_65 <- (threeC$age0-65)/10
  # threeC$age.last65 <- (threeC$age.last-65)/10
  # threeC$age.first65 <- (threeC$age.first-65)/10
  # threeC$SBP <- threeC$SBP/10
  # 
  # # arguments
  # formFixedScale = SBP ~ age.visit65
  # formRandomScale = ~ age.visit65
  # Location = 0
  # tau = runif(500, min = 0.2, max = 0.8)
  # formGroup = ~ ID
  # lambda = 0
  # data = threeC
  # n.chains = 3
  # n.iter = 600
  # n.burnin = 100
  # n.thin = 1
  # n.adapt = NULL
  # precision = 0.1
  # save_jagsUI = TRUE
  # parallel = FALSE

  # checks of arguments
  if(!(Location %in% c(0,1)))
    stop("'Location' is '0' or '1' to fix it to 0 or fit an intercept, respectively.")
  if(prod(tau<1)*prod(tau>0)!=1)
    stop("The 'tau' must be a scalar or a vector with value(s) in ]0;1[.")

  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),
                             all.vars(formFixedScale),
                             all.vars(formRandomScale)))]
  y <- data_long[all.vars(formFixedScale)][, 1]
  mfX <- stats::model.frame(formFixedScale, data = data_long)
  X2 <- stats::model.matrix(formFixedScale, mfX)
  mfZ <- stats::model.frame(formRandomScale, data = data_long)
  Z2 <- stats::model.matrix(formRandomScale, mfZ)
  ncX2 <- ncol(X2)
  ncZ2 <- ncol(Z2)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  if(length(tau)>1 && length(tau)!=I)
    stop("Length of `tau` must be 1 or egal to number of ID")
  if(length(tau)==1)
    tau <- rep(tau, I)

  # list of data jags
  jags.data <- list(y = y,
                    X2 = X2,
                    Z2 = Z2,
                    tau = tau,
                    ncX2 = ncX2,
                    ncZ2 = ncZ2,
                    I = I,
                    offset = offset)
  
  jags_code <- "model{
  # likelihood
  for (i in 1:I){
    # constants
    c1[i] <- (1-2*tau[i])/(tau[i]*(1-tau[i]))
    c2[i] <- 2/(tau[i]*(1-tau[i]))
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- beta + b[i] + c1[i]*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2[i]*va1[j])
    }#end of j loop
    # random effects
    b[i] ~ dmnorm(0, prec.Sigma_b)
    u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
  }#end of i loop
  # priors for parameters
  prec.Sigma_b ~ dgamma(priorA.b, priorB.b)
  covariance.b <- inverse(prec.Sigma_b)
  beta ~ dnorm(priorMean.beta, priorTau.beta)
  prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)
  covariance.u <- inverse(prec.Sigma_u[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
}"
  
  # manage scale part
  # regression on random effects u
  if(ncZ2==1){
    jags_code <- gsub("inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])",
                      "u[i] * Z2[j]",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])",
                      "u[i] ~ dnorm(mu_u, prec.Sigma_u)",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)",
                      "prec.Sigma_u ~ dgamma(priorA.u, priorB.u)",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("covariance.u <- inverse(prec.Sigma_u[, ])",
                      "covariance.u <- 1/prec.Sigma_u",
                      jags_code, fixed = TRUE)
    # initialization of chains
    initial.values <- list(u = rep(0, I),
                           prec.Sigma_u = precision)
    # complete jags data
    jags.data <- c(jags.data,
                   list(mu_u = 0,
                        priorA.u = precision,
                        priorB.u = precision)
                   )
    jags.data$ncZ2 <- NULL 
    jags.data$Z2 <- as.numeric(jags.data$Z2)
  }else{
    rplc <- paste(paste("u[i, ", 1:jags.data$ncZ2, "] * Z2[j, ", 1:jags.data$ncZ2, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])", rplc, jags_code, fixed = TRUE)
    # initialization of chains
    initial.values <- list(u = matrix(0, ncol = ncZ2, nrow = I),
                           prec.Sigma_u = diag(rep(precision, ncZ2)))
    # complete jags data
    jags.data <- c(jags.data,
                   list(mu_u = rep(0, ncZ2),
                        priorA.u = diag(rep(precision, ncZ2)),
                        priorB.u = ncZ2)
    )
  }
  if(ncX2==1){
    jags_code <- gsub("inprod(xi[1:ncX2], X2[j, 1:ncX2])", 
                      "xi * X2[j]", 
                      jags_code, fixed = TRUE)
    jags_code <- gsub("xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])",
                      "xi ~ dnorm(priorMean.xi, priorTau.xi)",
                      jags_code, fixed = TRUE)
    # initialization of chains
    initial.values <- c(initial.values,
                        list(xi = 0)
                        )
    # complete jags data
    jags.data <- c(jags.data,
                   list(priorMean.xi = 0,
                        priorTau.xi = precision))
    jags.data$ncX2 <- NULL
    jags.data$X2 <- as.numeric(jags.data$X2)
  }else{
    rplc <- paste(paste("xi[", 1:jags.data$ncX2, "] * X2[j, ", 1:jags.data$ncX2, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(xi[1:ncX2], X2[j, 1:ncX2])", rplc, jags_code, fixed = TRUE)
    # initialization of chains
    initial.values <- c(initial.values,
                        list(xi = rep(0, ncX2))
    )
    # complete jags data
    jags.data <- c(jags.data,
                   list(priorMean.xi = rep(0, ncX2),
                        priorTau.xi = diag(rep(precision, ncX2))))
  }
  if(Location==0){
    jags_code <- gsub("beta + b[i] +", 
                      "", 
                      jags_code, fixed = TRUE)
    jags_code <- gsub("b[i] ~ dmnorm(0, prec.Sigma_b)",
                      "",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("prec.Sigma_b ~ dgamma(priorA.b, priorB.b)
  covariance.b <- inverse(prec.Sigma_b)
  beta ~ dnorm(priorMean.beta, priorTau.beta)",
                      "",
                      jags_code, fixed = TRUE)
  }else{
    # initialization of chains
    initial.values <- c(initial.values,
                        list(b = rep(0, I),
                             beta = 0,
                             prec.Sigma_b = precision)
    )
    # complete jags data
    jags.data <- c(jags.data,
                   list(priorMean.beta = 0,
                        priorTau.beta = precision,
                        priorA.b = precision,
                        priorB.b = precision)
    )
  }
  
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
  parms_to_save <- c("xi", "u", "covariance.u")
  if(Location==1){
    parms_to_save <- c(parms_to_save, "beta", "b", "covariance.b")
  }


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
                          DIC = TRUE)

  #---- output
  out <- list(data = data)
  out$control <- list(formFixedScale = formFixedScale,
                      formRandomScale = formRandomScale,
                      formGroup = formGroup,
                      Location = Location,
                      tau = tau,
                      lambda = lambda,
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

  # b random effect output
  if(Location == 1){
    out$sims.list$b <- NULL
    out$random_effect_b <- list(postMeans = out_jags$mean$b,
                                postSd = out_jags$sd$b,
                                Rhat = out_jags$Rhat$b)
    names(out$random_effect_b$postMeans) <-
      names(out$random_effect_b$postSd) <-
      names(out$random_effect_b$Rhat) <- "(intercept)"
  }
  
  # u random effect output
  out$sims.list$u <- NULL
  out$random_effect_u <- list(postMeans = out_jags$mean$u,
                              postSd = out_jags$sd$u,
                              Rhat = out_jags$Rhat$u)
  if(ncZ2==1){
    names(out$random_effect_u$postMeans) <-
      names(out$random_effect_u$postSd) <-
      names(out$random_effect_u$Rhat) <- colnames(Z2)
  }else{
    colnames(out$random_effect_u$postMeans) <-
      colnames(out$random_effect_u$postSd) <-
      colnames(out$random_effect_u$Rhat) <- colnames(Z2)
  }

  # median : Posterior median of parameters
  out$median <- out_jags$q50
  if(Location == 1) out$median$b <- NULL
  out$median$u <- NULL

  # mean : Posterior mean of parameters
  out$mean <- out_jags$mean
  if(Location == 1) out$mean$b <- NULL
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
  if(Location == 1) out$StDev$b <- NULL
  out$StDev$u <- NULL

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  if(Location == 1) out$Rhat$b <- NULL
  out$Rhat$u <- NULL

  # names
  if(Location == 1){
    names(out$mean$beta) <-
      names(out$median$beta) <-
      names(out$modes$beta) <-
      names(out$StErr$beta) <-
      names(out$Rhat$beta) <-
      names(out$StDev$beta) <- "(intercept)"
    names(out$mean$covariance.b) <-
      names(out$median$covariance.b) <-
      names(out$modes$covariance.b) <-
      names(out$StErr$covariance.b) <-
      names(out$Rhat$covariance.b) <-
      names(out$StDev$covariance.b) <- "(intercept)"
  }

  names(out$mean$xi) <-
    names(out$median$xi) <-
    names(out$modes$xi) <-
    names(out$StErr$xi) <-
    names(out$Rhat$xi) <-
    names(out$StDev$xi) <- colnames(X2)

    if(ncZ2==1){
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

  # credible intervals

  ## beta parameters
  if(Location == 1){
    out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                          as.vector(t(out_jags$q97.5$beta)))
    rownames(out$CIs$beta) <- "(intercept)"
    colnames(out$CIs$beta) <- c("2.5%", "97.5%")
  }

  ## xi parameters
  out$CIs$xi <- cbind(as.vector(t(out_jags$q2.5$xi)),
                      as.vector(t(out_jags$q97.5$xi)))
  rownames(out$CIs$xi) <- colnames(X2)
  colnames(out$CIs$xi) <- c("2.5%", "97.5%")

  # only for diagonal elements of covariance matrix of random effects
  if(Location == 1){
    out$CIs$covariance.b <- cbind(as.vector(t(out_jags$q2.5$covariance.b)),
                                  as.vector(t(out_jags$q97.5$covariance.b)))
    rownames(out$CIs$covariance.b) <- "(intercept)"
    colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  }
  if(ncZ2==1){
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

  #---- End of the function defining the class and returning the output
  class(out) <- "BscALDreg"
  out

}
