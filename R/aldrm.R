#' \code{aldrm} fits asymmetric Laplace distribution regression mixed model for longitudinal data.
#'
#' Function using 'JAGS' sampler to estimate the individual response distribution assuming asymmetric Laplace
#' distribution.
#'
#' @param formFixed formula for fixed part of location regression with response variable
#' @param formRandom formula for random part of location regression without left formula part
#' @param formFixedScale formula for fixed part of scale regression without left formula part (log link)
#' @param formRandomScale formula for random part of scale regression without left formula part  (log link)
#' @param formFixedSkewness formula for fixed part of Skewness regression without left formula part (logit link)
#' @param formRandomSkewness formula for random part of Skewness regression without left formula part (logit link)
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param data long format dataset of observed variables
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param object_inits Object returns by 'lslqmm', 'lqmm' or 'lqm' functions of BeQuT package to initialize parameter chains, default is \code{NULL}
#' @param precision precision (e.g. inverse of variance) of prior distributions, default is 0.1
#' @param parallel see \code{jagsUI::jags()} function.
#'
#'
#' @return A \code{Baldrm} object is a list with the following elements:
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
#' #---- Fit ALD regression model with subject-specific random effect defining all distribution parameters
#' fit <- aldrm(formFixed = SBP ~ age.visit65,
#'              formRandom = ~ age.visit65,
#'              formFixedScale = ~ age.visit65,
#'              formRandomScale = ~ age.visit65,
#'              formFixedSkewness = ~ age.visit65,
#'              formRandomSkewness = ~ age.visit65,
#'              formGroup = ~ ID,
#'              data = threeC,
#'              n.iter = 3000,
#'              n.burnin = 1000)
#'
#' #---- Get the posterior statistics
#' fit$summary_table
#'
#' #---- Visualize the trace for parameters
#' MCMCvis::MCMCtrace(fit$jags$post,
#'                    params = c("beta", "xi", "alpha",
#'                               "covariance.b","covariance.u", "covariance.a"),
#'                    pdf = FALSE)
#'
#' #---- Plot of selected quantiles of the distribution over time
#' plot(fit, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), individuals = c(1,2,3,4), time.var = "age.visit65")
#'
#' }
#'
aldrm <- function(formFixed,
                  formRandom,
                  formFixedScale,
                  formRandomScale,
                  formFixedSkewness,
                  formRandomSkewness,
                  formGroup,
                  data,
                  n.chains = 3,
                  n.iter = 6000,
                  n.burnin = 1000,
                  n.thin = 1,
                  n.adapt = NULL,
                  object_inits = NULL,
                  precision = 0.1,
                  save_jagsUI = TRUE,
                  parallel = FALSE){


  #-- To do list
  # change the example
  # add the check on arguments
  # check when object_inits is given and add lqmm lslqmm

  # #---- Load data
  # data(threeC)
  # #data management
  # threeC$age.visit65 <- (threeC$age.visit-65)/10
  # threeC$age.final65 <- (threeC$age.final-65)/10
  # threeC$age0_65 <- (threeC$age0-65)/10
  # threeC$age.last65 <- (threeC$age.last-65)/10
  # threeC$age.first65 <- (threeC$age.first-65)/10
  # threeC$SBP <- threeC$SBP/10
  # # arguement
  # formFixed = SBP ~ age.visit65
  # formRandom = ~ age.visit65
  # formFixedScale = ~ age.visit65
  # formRandomScale = ~ age.visit65
  # formFixedSkewness = ~ age.visit65
  # formRandomSkewness = ~ age.visit65
  # formGroup = ~ ID
  # data = threeC
  # n.chains = 3
  # n.iter = 600
  # n.burnin = 100
  # n.thin = 1
  # n.adapt.step = 1000
  # n.adapt.max = 20000
  # n.adapt = NULL
  # object_inits = NULL
  # precision = 0.1
  # save_jagsUI = TRUE
  # parallel = FALSE

  # # #---- Load data
  # data(data_sim)
  # formFixed = y ~ visit + I(visit^2) + X1 + X2
  # formRandom = ~ visit + I(visit^2)
  # formFixedScale = ~ visit + X1
  # formRandomScale = ~ visit
  # formFixedSkewness = ~ -1 + X1 + X2
  # formRandomSkewness = ~ visit
  # formGroup = ~ id
  # data = data_sim$data_long
  # n.chains = 3
  # n.iter = 1000
  # n.burnin = 500
  # n.thin = 5
  # n.adapt.step = 1000
  # n.adapt.max = 20000
  # n.adapt = NULL
  # object_inits = NULL
  # precision = 0.1
  # save_jagsUI = TRUE
  # parallel = FALSE

  # checks of arguments
  if(!is.null(object_inits) && !inherits(object_inits, "Blslqmm"))
    stop("The 'object_inits' object must be an \"Blqm\" objects.")
  if(!is.null(object_inits) && formFixed!=object_inits$control$formFixed)
    stop("The 'formFixed' formula does not match the one from 'object_inits'.")

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
  if(!is.null(object_inits) && inherits(object_inits, "Blslqmm")){
    # prior beta parameters
    priorMean.beta <- as.vector(object_inits$mean$beta)
    priorTau.beta <- if(ncX1==1) c(precision) else diag((1/object_inits$StDev$beta^2)*precision)
    # prior xi parameters
    priorMean.xi <- as.vector(object_inits$mean$xi)
    priorTau.xi <- if(ncX2==1) c(precision) else diag((1/object_inits$StDev$xi^2)*precision)

  }else{
    # prior beta parameters
    priorMean.beta <- rep(0, ncX1)
    priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
    # prior xi parameters
    priorMean.xi <- if(ncX2==1) c(0) else rep(0, ncX2)
    priorTau.xi <- if(ncX2==1) c(precision) else diag(rep(precision, length(priorMean.xi)))
  }
  # initialization of chains
  initial.values <- list(b = if(is.null(object_inits)) matrix(0, ncol = ncZ1, nrow = I) else object_inits$random_effect$postMeans,
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
    if(is.null(object_inits)){
      initial.values$prec.Sigma_b <- diag(precision, ncol(Z1))
      initial.values$prec.Sigma_u <- diag(precision, ncol(Z2))
    }else{
      initial.values$prec.Sigma_b <- solve(object_inits$mean$covariance.b)
      initial.values$prec.Sigma_u <- solve(object_inits$mean$covariance.u)
    }


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

    #-----------------------------------------------------------------------------
    #---- Model definition step

    ## Large writting of aldrm
      jags_code <- "model{
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      # ald parameters
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])
      sigma[j] <- exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2]))
      logit(tau[j]) <- inprod(alpha[1:ncX3], X3[j, 1:ncX3]) + inprod(a[i, 1:ncZ3], Z3[j, 1:ncZ3])
      # mixture representation
      mu_N[j]  <- mu[j] + c1[j]*va1[j]
      prec[j] <- 1/(sigma[j]*c2[j]*va1[j])
      c1[j] <- (1-2*tau[j])/(tau[j]*(1-tau[j]))
      c2[j] <- 2/(tau[j]*(1-tau[j]))
      # distributions
      y[j] ~ dnorm(mu_N[j], prec[j])
      va1[j] ~ dexp(1/sigma[j])
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


    # #loglik_cond
    # loglik_cond[j] <- logdensity.norm(y[j], mu[j], prec[j])
    # #loglik_marg
    # resid[j] <- y[j] - mu[j]
    # rho[j] <- resid[j]*(tau[j] - 1 + step(resid[j]))
    # loglik_marg[j] <- log(tau[j]) + log(1 - tau[j]) - log(sigma[j]) - rho[j]/sigma[j]

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

  # model specification for skewness part
  if(ncZ3==1){
    jags.data <- c(jags.data,
                   list(priorA.a = precision,
                        priorB.a = precision,
                        mu_a = 0
                   ))
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
  }else{
      jags.data <- c(jags.data,
                     list(Z3 = Z3,
                          ncZ3 = ncZ3,
                          priorA.a = diag(rep(1, ncZ3)),
                          priorB.a = ncZ3,
                          mu_a = rep(0, ncZ3)
                     ))
      # regression on random effects a
      rplc <- paste(paste("a[i, ", 1:jags.data$ncZ3, "] * Z3[j, ", 1:jags.data$ncZ3, "]", sep = ""), collapse = " + ")
      jags_code <- gsub("inprod(a[i, 1:ncZ3], Z3[j, 1:ncZ3])", rplc, jags_code, fixed = TRUE)
  }
  if(ncX3==1){
    jags.data <- c(jags.data,
                        list(priorTau.alpha = 1,
                             priorMean.alpha = 0)
                   )

  }else{
    jags.data <- c(jags.data,
                   list(X3 = X3,
                        ncX3 = ncX3,
                        priorMean.alpha = rep(0,ncX3),
                        priorTau.alpha = as.matrix(priorTau.alpha)
                   ))
    # regression from X3 with alpha
    rplc <- paste(paste("alpha[", 1:jags.data$ncX3, "] * X3[j, ", 1:jags.data$ncX3, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(alpha[1:ncX3], X3[j, 1:ncX3])", rplc, jags_code, fixed = TRUE)
  }
  # if(WAIC){
  #   jags_code <- gsub("#loglik_cond",
  #                     "loglik_cond[j] <- logdensity.norm(y[j], mu[j], prec[j])",
  #                     jags_code, fixed = TRUE)
  #   jags_code <- gsub("#loglik_marg",
  #                     "resid[j] <- y[j] - mu[j]
  #                      rho[j] <- resid[j]*(tau[j] - 1 + step(resid[j]))
  #                      loglik_marg[j] <- log(tau[j]) + log(1 - tau[j]) - log(sigma[j]) - rho[j]/sigma[j]",
  #                     jags_code, fixed = TRUE)
  #   # where 1_{v<0} = 1 - step(v)
  # }


  # initialization values
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
  parms_to_save <- c("beta", "xi", "alpha",
                     "covariance.b", "covariance.u", "covariance.a",
                     "b", "u", "a")

  if(ncX3==1 && ncZ3==1){
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
                      formFixedScale = formFixedScale,
                      formRandomScale = formRandomScale,
                      formFixedSkewness = formFixedSkewness,
                      formRandomSkewness = formRandomSkewness,
                      formGroup = formGroup,
                      call_function = "aldrm",
                      n.chains = n.chains,
                      parallel = parallel,
                      n.iter = n.iter,
                      n.burnin = n.burnin,
                      n.thin = n.thin,
                      I = I)

  out$jags <- list(jags_code = jags_code,
                   jags.data = jags.data,
                   out_jags = out_jags)


  #- summary outputs
  parms_to_save2 <- c("beta", "xi", "alpha",
                      "covariance.b", "covariance.u", "covariance.a")

  out$summary_table <- MCMCvis::MCMCsummary(out$jags$out_jags$samples,
                                            params = parms_to_save2)


  #---- End of the function defining the class and returning the output
  class(out) <- c("Baldrm", "Bald", "Bdr", "Bmm")
  out

}
