#' \code{ALDreg} fits density regression for longitudinal data assuming asymmetric Laplace distribution
#'
#' Function using 'JAGS' software to estimate the individual response distribution assuming asymmetric Laplace
#' distribution.
#'
#' @param formFixed formula for fixed part of location regression with response variable
#' @param formRandom formula for random part of location regression without left formula part
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param formFixedScale formula for fixed part of scale regression without left formula part (log link)
#' @param formRandomScale formula for random part of scale regression without left formula part  (log link)
#' @param formFixedSkewness formula for fixed part of Skewness regression without left formula part (logit link)
#' @param formRandomSkewness formula for random part of Skewness regression without left formula part (logit link)
#' @param Skewness_link Considered link between the linear predictor and the skewness parameter. By default "logit", but considering only group-specific constant, the beta distribution can use as prior by specifying "beta".
#' @param data long format dataset of observed variables
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param object_lslqmm Object returns by 'lslqmm' function of BeQuT package to initialize parameter chains, default is \code{NULL}
#' @param precision precision (e.g. inverse of variance) of prior distributions, default is 0.1
#' @param save_jagsUI If \code{TRUE} (by default), the output of \code{jagsUI} package is returned by the function. Warning, if \code{TRUE}, the output can be large.
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{BALDreg} object is a list with the following elements:
#'  \describe{
#'   \item{\code{mean}}{list of posterior mean for each parameter}
#'   \item{\code{median}}{list of posterior median for each parameter}
#'   \item{\code{modes}}{list of posterior mode for each parameter}
#'   \item{\code{StErr}}{list of standard error for each parameter}
#'   \item{\code{StDev}}{list of standard deviation for each parameter}
#'   \item{\code{Rhat}}{list of Rhat convergence criteria for each parameter}
#'   \item{\code{ICs}}{list of the credibility interval at 0.95 for each parameters excepted for covariance parameters in covariance matrix of random effects. Otherwise, use save_jagsUI=TRUE to have the associated quantiles.}
#'   \item{\code{data}}{data included in argument}
#'   \item{\code{sims.list}}{list of the MCMC chains of the parameters and random effects}
#'   \item{\code{control}}{list of arguments giving details about the estimation}
#'   \item{\code{random_effect_b}}{list including posterior mean, posterior standard deviation and Rhat criteria of subject-specific location random effects}
#'   \item{\code{random_effect_u}}{list including posterior mean, posterior standard deviation and Rhat criteria of subject-specific scale random effects}
#'   \item{\code{random_effect_a}}{list including posterior mean, posterior standard deviation and Rhat criteria of subject-specific skewness random effects}
#'   \item{\code{random_effect_tau}}{only if no covariates are considered in skewness part, list including posterior mean, posterior standard deviation and Rhat criteria of subject-specific skewness random effect}
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
#' #---- Fit ALD regression model with subject-specific random effect for skewness part
#' model_beta <- ALDreg(formFixed = SBP ~ age.visit65,
#'                         formRandom = ~ age.visit65,
#'                         formFixedScale = ~ age.visit65,
#'                         formRandomScale = ~ age.visit65,
#'                         formGroup = ~ ID,
#'                         Skewness_link = "beta",
#'                         data = threeC,
#'                         n.iter = 3000,
#'                         n.burnin = 1000)
#'
#' #---- Fit ALD regression model with subject-specific random effect for skewness part
#' model_intercept <- ALDreg(formFixed = SBP ~ age.visit65,
#'                         formRandom = ~ age.visit65,
#'                         formFixedScale = ~ age.visit65,
#'                         formRandomScale = ~ age.visit65,
#'                         formGroup = ~ ID,
#'                         Skewness_link = "logit",
#'                         data = threeC,
#'                         n.iter = 3000,
#'                         n.burnin = 1000)
#'
#' #---- Get the posterior means
#' model_intercept$mean
#'
#' #---- Get the posterior prediction of subject-specific skewness part
#' model_intercept$out$random_effect_tau$postMeans
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(model_intercept$out_jagsUI, parameters = "beta")
#'
#' #---- Summary of output
#' summary(model_intercept)
#'
#' #---- Fit full ALD regression model with subject-specific random effect defining all parts
#' model_full <- ALDreg(formFixed = SBP ~ age.visit65,
#'                      formRandom = ~ age.visit65,
#'                      formFixedScale = ~ age.visit65,
#'                      formRandomScale = ~ age.visit65,
#'                      formFixedSkewness = ~ age.visit65,
#'                      formRandomSkewness = ~ age.visit65,
#'                      formGroup = ~ ID,
#'                      data = threeC,
#'                      n.iter = 3000,
#'                      n.burnin = 1000)
#' }
#'
ALDreg <- function(formFixed,
                   formRandom,
                   formGroup,
                   formFixedScale,
                   formRandomScale,
                   formFixedSkewness = NULL,
                   formRandomSkewness = NULL,
                   Skewness_link = "logit",
                   data,
                   n.chains = 3,
                   n.iter = 6000,
                   n.burnin = 1000,
                   n.thin = 1,
                   n.adapt = NULL,
                   object_lslqmm = NULL,
                   precision = 0.1,
                   save_jagsUI = TRUE,
                   parallel = FALSE){


  #-- To do list
  # change the example
  # add the check on arguments
  # check when object_B is given and add lqmm lslqmm

  # #---- Load data
  # data(threeC)
  #
  # #---- data management
  # threeC$age.visit65 <- (threeC$age.visit-65)/10
  # threeC$age.final65 <- (threeC$age.final-65)/10
  # threeC$age0_65 <- (threeC$age0-65)/10
  # threeC$age.last65 <- (threeC$age.last-65)/10
  # threeC$age.first65 <- (threeC$age.first-65)/10
  # threeC$SBP <- threeC$SBP/10
  #
  #
  #
  # formFixed = SBP ~ age.visit65
  # formRandom = ~ age.visit65
  # formFixedScale = ~ age.visit65
  # formRandomScale = ~ age.visit65
  # formFixedSkewness = ~ age.visit65
  # formRandomSkewness = ~ age.visit65
  # formGroup = ~ ID
  # Skewness_link = "beta"
  # data = threeC
  # n.chains = 3
  # n.iter = 600
  # n.burnin = 100
  # n.thin = 1
  # n.adapt = NULL
  # object_lslqmm = NULL
  # precision = 0.1
  # save_jagsUI = TRUE
  # parallel = FALSE

  # checks of arguments
  if(!is.null(object_lslqmm) && !inherits(object_lslqmm, "Blslqmm"))
    stop("The 'object_lslqmm' obejct must be an \"Blslqmm\" objects.")
  if(!is.null(object_lslqmm) && formFixed!=object_lslqmm$control$formFixed)
    stop("The 'formFixed' formula does not match the one from 'object_lslqmm'.")
  if(Skewness_link != "logit")
    formFixedSkewness <- formRandomSkewness <- NULL

  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),
                             all.vars(formFixed),
                             all.vars(formRandom),
                             all.vars(formFixedScale),
                             all.vars(formRandomScale),
                             all.vars(formFixedSkewness),
                             all.vars(formRandomSkewness)))]
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
  if(!is.null(object_lslqmm) && inherits(object_lslqmm, "Blslqmm")){
    # prior beta parameters
    priorMean.beta <- as.vector(object_lslqmm$mean$beta)
    priorTau.beta <- if(ncX1==1) c(precision) else diag((1/object_lslqmm$StDev$beta^2)*precision)
    # prior xi parameters
    priorMean.xi <- as.vector(object_lslqmm$mean$xi)
    priorTau.xi <- if(ncX2==1) c(precision) else diag((1/object_lslqmm$StDev$xi^2)*precision)

  }else{
    # prior beta parameters
    priorMean.beta <- rep(0, ncX1)
    priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
    # prior xi parameters
    priorMean.xi <- if(ncX2==1) c(0) else rep(0, ncX2)
    priorTau.xi <- if(ncX2==1) c(precision) else diag(rep(precision, length(priorMean.xi)))
  }
  # initialization of chains
  initial.values <- list(b = if(is.null(object_lslqmm)) matrix(0, ncol = ncZ1, nrow = I) else object_lslqmm$random_effect$postMeans,
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
    if(is.null(object_lslqmm)){
      initial.values$prec.Sigma_b <- diag(precision, ncol(Z1))
      initial.values$prec.Sigma_u <- diag(precision, ncol(Z2))
    }else{
      initial.values$prec.Sigma_b <- solve(object_lslqmm$mean$covariance.b)
      initial.values$prec.Sigma_u <- solve(object_lslqmm$mean$covariance.u)
    }

    # define the appropriate BUGS model
    ## tau is a group-specific parameter
    if(Skewness_link == "beta"){
      jags_code <- "model{
  # likelihood
  for (i in 1:I){
    # constants
    c1[i] <- (1-2*tau[i])/(tau[i]*(1-tau[i]))
    c2[i] <- 2/(tau[i]*(1-tau[i]))
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1[i]*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2[i]*va1[j])
    }#end of j loop
    # random effects
    b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])
    u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
    tau[i] ~ dbeta(1,1)
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

    if(Skewness_link == "logit"){
      if(is.null(formFixedSkewness) && is.null(formRandomSkewness))
        formFixedSkewness <- formRandomSkewness <- ~ 1
      mfX <- stats::model.frame(formFixedSkewness, data = data_long)
      X3 <- stats::model.matrix(formFixedSkewness, mfX)
      mfZ <- stats::model.frame(formRandomSkewness, data = data_long)
      Z3 <- stats::model.matrix(formRandomSkewness, mfZ)
      ncX3 <- ncol(X3)
      ncZ3 <- ncol(Z3)
      # Initialization of prior parameter
      # prior alpha parameters
      priorMean.alpha <- if(ncX3==1) c(0) else rep(0, ncX3)
      priorTau.alpha <- if(ncX3==1) c(1) else diag(rep(1, length(priorMean.alpha)))
      # initialization of chains
      initial.values <- c(initial.values,
                          list(a = if(ncZ3==1) rep(0, I) else matrix(0, ncol = ncZ3, nrow = I),
                               alpha = priorMean.alpha
                               )
                          )

      # define the jags model
      if(ncX3==1 && ncZ3==1){
        jags_code <- "model{
  # likelihood
  for (i in 1:I){
    # constants
    c1[i] <- (1-2*tau[i])/(tau[i]*(1-tau[i]))
    c2[i] <- 2/(tau[i]*(1-tau[i]))
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1[i]*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2[i]*va1[j])
    }#end of j loop
    # random effects
    b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])
    u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
    logit(tau[i]) <- alpha + a[i]
    a[i] ~ dnorm(mu_a, prec.Sigma_a)
  }#end of i loop
  # priors for parameters
  prec.Sigma_b[1:ncZ1, 1:ncZ1] ~ dwish(priorA.b[, ], priorB.b)
  covariance.b <- inverse(prec.Sigma_b[, ])
  prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)
  covariance.u <- inverse(prec.Sigma_u[, ])
  prec.Sigma_a ~ dgamma(priorA.a, priorB.a)
  covariance.a <- 1/prec.Sigma_a
  beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
  alpha ~ dnorm(priorMean.alpha, priorTau.alpha)
}"
        jags.data <- c(jags.data,
                       list(priorTau.alpha = 1,
                            priorMean.alpha = 0,
                            priorA.a = 1,
                            priorB.a = 1,
                            mu_a = 0
                            )
                       )
      }
#       if(ncX3>1 && ncZ3==1){
#         jags_code <- "model{
#   # likelihood
#   for (i in 1:I){
#     # longitudinal part
#     for(j in offset[i]:(offset[i+1]-1)){
#       y[j] ~ dnorm(mu[j], prec[j])
#       mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1[j]*va1[j]
#       sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
#       va1[j] ~ dexp(1/sigma[j])
#       prec[j] <- 1/(sigma[j]*c2[j]*va1[j])
#       # constants
#       c1[j] <- (1-2*tau[j])/(tau[j]*(1-tau[j]))
#       c2[j] <- 2/(tau[j]*(1-tau[j]))
#       logit(tau[j]) <- inprod(alpha[1:ncX3], X3[j, 1:ncX3]) + a[i]
#     }#end of j loop
#     # random effects
#     b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])
#     u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
#     a[i] ~ dnorm(mu_a, priorA.a)
#   }#end of i loop
#   # priors for parameters
#   prec.Sigma_b[1:ncZ1, 1:ncZ1] ~ dwish(priorA.b[, ], priorB.b)
#   covariance.b <- inverse(prec.Sigma_b[, ])
#   prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)
#   covariance.u <- inverse(prec.Sigma_u[, ])
#   beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
#   xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
#   alpha ~ dnorm(priorMean.alpha, priorTau.alpha)
# }"
#       }
      # if(ncX3>1 && ncZ3>1){
      if(ncX3>1){
        jags_code <- "model{
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1]) + c1[j]*va1[j]
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      va1[j] ~ dexp(1/sigma[j])
      prec[j] <- 1/(sigma[j]*c2[j]*va1[j])
      # constants
      c1[j] <- (1-2*tau[j])/(tau[j]*(1-tau[j]))
      c2[j] <- 2/(tau[j]*(1-tau[j]))
      logit(tau[j]) <- inprod(alpha[1:ncX3], X3[j, 1:ncX3]) + inprod(a[i, 1:ncZ3], Z3[j, 1:ncZ3])
    }#end of j loop
    # random effects
    b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])
    u[i, 1:ncZ2] ~ dmnorm(mu_u[], prec.Sigma_u[, ])
    a[i, 1:ncZ3] ~ dmnorm(mu_a[], prec.Sigma_a[, ])
  }#end of i loop
  # priors for parameters
  prec.Sigma_b[1:ncZ1, 1:ncZ1] ~ dwish(priorA.b[, ], priorB.b)
  covariance.b <- inverse(prec.Sigma_b[, ])
  prec.Sigma_u[1:ncZ2, 1:ncZ2] ~ dwish(priorA.u[, ], priorB.u)
  covariance.u <- inverse(prec.Sigma_u[, ])
  prec.Sigma_a[1:ncZ3, 1:ncZ3] ~ dwish(priorA.a[, ], priorB.a)
  covariance.a <- inverse(prec.Sigma_a[, ])
  beta[1:ncX1] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  xi[1:ncX2] ~ dmnorm(priorMean.xi[], priorTau.xi[, ])
  alpha[1:ncX3] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
}"

        jags.data <- c(jags.data,
                       list(X3 = X3,
                            ncX3 = ncX3,
                            priorMean.alpha = rep(0,ncX3),
                            priorTau.alpha = as.matrix(priorTau.alpha)
                       ))
        }
      }


  }else{
    stop("ALDreg is implemented only when both location and scale design vectors are a length >1.")
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
  # regression on random effects u
  rplc <- paste(paste("u[i, ", 1:jags.data$ncZ2, "] * Z2[j, ", 1:jags.data$ncZ2, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])", rplc, jags_code, fixed = TRUE)
  if(Skewness_link == "logit" && ncX3>1 && ncZ3>1){
    jags.data <- c(jags.data,
                   list(Z3 = Z3,
                        ncZ3 = ncZ3,
                        priorA.a = diag(rep(1, ncZ3)),
                        priorB.a = ncZ3,
                        mu_a = rep(0, ncZ3)
                   ))
    # regression from X3 with alpha
    rplc <- paste(paste("alpha[", 1:jags.data$ncX3, "] * X3[j, ", 1:jags.data$ncX3, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(alpha[1:ncX3], X3[j, 1:ncX3])", rplc, jags_code, fixed = TRUE)
    # regression on random effects a
    rplc <- paste(paste("a[i, ", 1:jags.data$ncZ3, "] * Z3[j, ", 1:jags.data$ncZ3, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(a[i, 1:ncZ3], Z3[j, 1:ncZ3])", rplc, jags_code, fixed = TRUE)
  }
  if(Skewness_link == "logit" && ncX3>1 && ncZ3==1){
    # regression from X3 with alpha
    rplc <- paste(paste("alpha[", 1:jags.data$ncX3, "] * X3[j, ", 1:jags.data$ncX3, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(alpha[1:ncX2], X3[j, 1:ncX3])", rplc, jags_code, fixed = TRUE)
    # regression on random effects a
    jags_code <- gsub("inprod(a[i, 1:ncZ3], Z3[j, 1:ncZ3])",
                      "a[i]",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("a[i, 1:ncZ3] ~ dmnorm(mu_a[], prec.Sigma_a[, ])",
                      "a[i] ~ dnorm(mu_a, prec.Sigma_a)",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("prec.Sigma_a[1:ncZ3, 1:ncZ3] ~ dwish(priorA.a[, ], priorB.a)",
                      "prec.Sigma_a ~ dgamma(priorA.a, priorB.a)",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("covariance.a <- inverse(prec.Sigma_a[, ])",
                      "covariance.a <- 1/prec.Sigma_a",
                      jags_code, fixed = TRUE)
    jags.data <- c(jags.data,
                   list(priorA.a = 1,
                        mu_a = 0
                   ))
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
  parms_to_save <- c("beta", "xi", "b", "covariance.b","u", "covariance.u")
  if(Skewness_link == "beta"){
    parms_to_save <- c(parms_to_save, "tau")
  }
  if(Skewness_link == "logit"){
    parms_to_save <- c(parms_to_save, "alpha", "a", "covariance.a")
  }
  if(Skewness_link == "logit" && ncX3==1 && ncZ3==1){
    parms_to_save <- c(parms_to_save, "tau")
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
  out$control <- list(formFixed = formFixed,
                      formRandom = formRandom,
                      formGroup = formGroup,
                      formFixedScale = formFixedScale,
                      formRandomScale = formRandomScale,
                      formFixedSkewness = formFixedSkewness,
                      formRandomSkewness = formRandomSkewness,
                      Skewness_link = Skewness_link,
                      call_function = "ALDreg",
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

  # tau or a random effect output
  if(Skewness_link == "logit"){
    out$sims.list$a <- NULL
    out$random_effect_a <- list(postMeans = out_jags$mean$a,
                                postSd = out_jags$sd$a,
                                Rhat = out_jags$Rhat$a)
    names(out$random_effect_a$postMeans) <-
             names(out$random_effect_a$postSd) <-
             names(out$random_effect_a$Rhat) <- if(ncZ3==1) "intercept" else colnames(Z3)

  }
  if(Skewness_link == "beta" || (Skewness_link == "logit" && ncX3==1 && ncZ3==1)){
    out$sims.list$tau <- NULL
    out$random_effect_tau <- list(postMeans = out_jags$mean$tau,
                                  postSd = out_jags$sd$tau,
                                  Rhat = out_jags$Rhat$tau)
  }

  # median : Posterior median of parameters
  out$median <- out_jags$q50
  out$median$b <- NULL
  out$median$u <- NULL
  if(Skewness_link == "beta" || (Skewness_link == "logit" && ncX3==1 && ncZ3==1))
    out$median$tau <- NULL
  if(Skewness_link == "logit")
    out$median$a <- NULL

  # mean : Posterior mean of parameters
  out$mean <- out_jags$mean
  out$mean$b <- NULL
  out$mean$u <- NULL
  if(Skewness_link == "beta" || (Skewness_link == "logit" && ncX3==1 && ncZ3==1))
    out$mean$tau <- NULL
  if(Skewness_link == "logit")
    out$mean$a <- NULL

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
  if(Skewness_link == "beta" || (Skewness_link == "logit" && ncX3==1 && ncZ3==1))
    out$StDev$tau <- NULL
  if(Skewness_link == "logit")
    out$StDev$a <- NULL

  # Rhat : Gelman & Rubin diagnostic
  out$Rhat <- out_jags$Rhat
  out$Rhat$b <- NULL
  out$Rhat$u <- NULL
  if(Skewness_link == "beta" || (Skewness_link == "logit" && ncX3==1 && ncZ3==1))
    out$Rhat$tau <- NULL
  if(Skewness_link == "logit")
    out$Rhat$a <- NULL

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

  if(Skewness_link == "logit")
    names(out$mean$alpha) <-
    names(out$median$alpha) <-
    names(out$modes$alpha) <-
    names(out$StErr$alpha) <-
    names(out$Rhat$alpha) <-
    names(out$StDev$alpha) <- colnames(X3)

  if(ncZ1==1){
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

  if(Skewness_link == "logit"){
    if(ncZ3==1){
      names(out$mean$covariance.a) <-
        names(out$median$covariance.a) <-
        names(out$modes$covariance.a) <-
        names(out$StErr$covariance.a) <-
        names(out$Rhat$covariance.a) <-
        names(out$StDev$covariance.a) <- colnames(Z3)
    }else{
      colnames(out$mean$covariance.a) <-
        rownames(out$mean$covariance.a) <-
        colnames(out$median$covariance.a) <-
        rownames(out$median$covariance.a) <-
        colnames(out$modes$covariance.a) <-
        rownames(out$modes$covariance.a) <-
        colnames(out$Rhat$covariance.a) <-
        rownames(out$Rhat$covariance.a) <-
        colnames(out$StErr$covariance.a) <-
        rownames(out$StErr$covariance.a) <-
        colnames(out$StDev$covariance.a) <-
        rownames(out$StDev$covariance.a) <- colnames(Z3)
    }
  }

  # credible intervals

  ## beta parameters
  out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
                        as.vector(t(out_jags$q97.5$beta)))
  rownames(out$CIs$beta) <- colnames(X1)
  colnames(out$CIs$beta) <- c("2.5%", "97.5%")

  ## xi parameters
  out$CIs$xi <- cbind(as.vector(t(out_jags$q2.5$xi)),
                      as.vector(t(out_jags$q97.5$xi)))
  rownames(out$CIs$xi) <- colnames(X2)
  colnames(out$CIs$xi) <- c("2.5%", "97.5%")

  # only for diagonal elements of covariance matrix of random effects
  if(ncZ1==1){
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
  if(Skewness_link == "logit"){
    ## alpha parameters
    out$CIs$alpha <- cbind(as.vector(t(out_jags$q2.5$alpha)),
                           as.vector(t(out_jags$q97.5$alpha)))
    rownames(out$CIs$alpha) <- colnames(X3)
    colnames(out$CIs$alpha) <- c("2.5%", "97.5%")
    if(ncZ3==1){
      out$CIs$covariance.a <- cbind(as.vector(t(out_jags$q2.5$covariance.a)),
                                    as.vector(t(out_jags$q97.5$covariance.a)))
      rownames(out$CIs$covariance.a) <- colnames(Z3)
      colnames(out$CIs$covariance.a) <- c("2.5%", "97.5%")
    }else{
      out$CIs$covariance.a <- cbind(as.vector(diag(out_jags$q2.5$covariance.a)),
                                    as.vector(diag(out_jags$q97.5$covariance.a)))
      rownames(out$CIs$covariance.a) <- colnames(Z3)
      colnames(out$CIs$covariance.a) <- c("2.5%", "97.5%")
    }
  }

  # save jags output if requires
  if(save_jagsUI)
    out$out_jagsUI <- out_jags

  #---- End of the function defining the class and returning the output
  class(out) <- "BALDreg"
  out

}
