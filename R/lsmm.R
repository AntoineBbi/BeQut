#' \code{lsmm} fits location-scale mixed model
#'
#' Function using 'JAGS' software to estimate a location-scale mixed model (lsmm) assuming Gaussian
#' distribution for residual error.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without left formula part
#' @param formFixedScale formula for fixed part of scale regression part without left formula part
#' @param formRandomScale formula for random part of scale regression part without left formula part
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param data long format dataset of observed variables
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param precision precision (e.g. inverse of variance) of prior distributions, default is 0.1
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{Blsmm} object is a list with the following elements:
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
#' @references Courcoul L, Tzourio C, Woodward M, Barbieri A, Jacqmin-Gadda H (2024).
#' \emph{	A location-scale joint model for studying the link between the time-dependent subject-specific variability of blood pressure and competing events}.
#' arXiv, arXiv:2306.16785. doi: https://doi.org/10.48550/arXiv.2306.16785.
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
#' fit_lsmm <- lsmm(formFixed = SBP ~ age.visit65,
#'                         formRandom = ~ age.visit65,
#'                         formFixedScale = ~ age.visit65,
#'                         formRandomScale = ~ age.visit65,
#'                         formGroup = ~ ID,
#'                         data = threeC,
#'                         n.iter = 6000,
#'                         n.burnin = 1000)
#'
#' #---- Get the posterior statistics
#' fit_lsmm$summary_table
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(fit_lsmm$out_jagsUI, parameters = "beta")
#'
#' #---- Summary of output
#' summary(fit_lsmm)
#' }
#'
lsmm <- function(formFixed,
                   formRandom,
                   formGroup,
                   formFixedScale,
                   formRandomScale,
                   data,
                   n.chains = 3,
                   n.iter = 10000,
                   n.burnin = 5000,
                   n.thin = 1,
                   n.adapt = NULL,
                   precision = 0.1,
                   parallel = FALSE){


  #-- To do list
  # change the example
  # add the check on arguments

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


  priorMean.beta <- if(ncX1==1) c(0) else rep(0, ncX1)
  priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
  # prior xi parameters
  priorMean.xi <- if(ncX2==1) c(0) else rep(0, ncX2)
  priorTau.xi <- if(ncX2==1) c(precision) else diag(rep(precision, length(priorMean.xi)))

  # initialization of chains
  initial.values <- list(b = matrix(0, ncol = ncZ1, nrow = I),
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


  # define the appropriate BUGS model
  jags_code <- "model{
  # likelihood
  for (i in 1:I){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      y[j] ~ dnorm(mu[j], prec[j])
      mu[j] <- inprod(beta[1:ncX1], X1[j, 1:ncX1]) + inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])
      prec[j] <- pow(exp(inprod(xi[1:ncX2], X2[j, 1:ncX2]) + inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])), -2)
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

  # regression from X1 with beta
  rplc <- paste(paste("beta[", 1:jags.data$ncX1, "] * X1[j, ", 1:jags.data$ncX1, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX1], X1[j, 1:ncX1])", rplc, jags_code, fixed = TRUE)
  # regression on random effects b
  if(ncZ1>1){
  rplc <- paste(paste("b[i, ", 1:jags.data$ncZ1, "] * Z1[j, ", 1:jags.data$ncZ1, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])", rplc, jags_code, fixed = TRUE)
  # data jags
  jags.data <- c(jags.data,
                 list(priorA.b = diag(rep(precision, ncZ1)),
                      priorB.b = ncZ1,
                      mu_b = rep(0, ncZ1)
                 )
  )
  # update initialisation values
  initial.values$prec.Sigma_b <- diag(precision, ncol(Z1))
  }else{ # ncZ1==1
    # data jags
    jags.data <- c(jags.data,
                   list(priorA.b = precision,
                        priorB.b = precision,
                        mu_b = 0
                        )
                   )
    jags.data$ncZ1 <- NULL
    # regression on random effects b
    jags_code <- gsub("inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])",
                      "b[i]",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("b[i, 1:ncZ1] ~ dmnorm(mu_b[], prec.Sigma_b[, ])",
                      "b[i] ~ dnorm(mu_b, prec.Sigma_b)",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("prec.Sigma_b[1:ncZ1, 1:ncZ1] ~ dwish(priorA.b[, ], priorB.b)",
                      "prec.Sigma_b ~ dgamma(priorA.b, priorB.b)",
                      jags_code, fixed = TRUE)
    jags_code <- gsub("covariance.b <- inverse(prec.Sigma_b[, ])",
                      "covariance.b <- 1/prec.Sigma_b",
                      jags_code, fixed = TRUE)

    # update initialisation values
    initial.values$prec.Sigma_b <-precision
  }

  # regression from X2 with xi
  rplc <- paste(paste("xi[", 1:jags.data$ncX2, "] * X2[j, ", 1:jags.data$ncX2, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(xi[1:ncX2], X2[j, 1:ncX2])", rplc, jags_code, fixed = TRUE)
  # regression on random effects u
  if(ncZ2>1){
    rplc <- paste(paste("u[i, ", 1:jags.data$ncZ2, "] * Z2[j, ", 1:jags.data$ncZ2, "]", sep = ""), collapse = " + ")
    jags_code <- gsub("inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])", rplc, jags_code, fixed = TRUE)
    # data jags
    jags.data <- c(jags.data,
                   list(priorA.u = diag(rep(precision, ncZ2)),
                        priorB.u = ncZ2,
                        mu_u = rep(0, ncZ2)
                   )
    )
    # update initialisation values
    initial.values$prec.Sigma_u <- diag(precision, ncol(Z2))
  }else{ # ncZ2==1
    # data jags
    jags.data <- c(jags.data,
                   list(priorA.u = precision,
                        priorB.u = precision,
                        mu_u = 0
                   )
    )
    jags.data$ncZ2 <- NULL
    # regression on random effects b
    jags_code <- gsub("inprod(u[i, 1:ncZ2], Z2[j, 1:ncZ2])",
                      "u[i]",
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
    # update initialisation values
    initial.values$prec.Sigma_u <-precision
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

  #---- use JAGS sampler via jagsUI
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = textConnection(jags_code), # cat(jags_code)
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
                      call_function = "lsmm",
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


  #- summary outputs
  parms_to_save2 <- c("beta", "xi",
                      "covariance.b", "covariance.u")
  out$summary_table <- MCMCvis::MCMCsummary(out$jags$out_jags$samples,
                                            params = parms_to_save2)



  # #- other outputs
  #
  # # sims.list output
  # out$sims.list <- out_jags$sims.list
  # out$sims.list$b <- NULL
  # out$sims.list$u <- NULL
  #
  # # b random effect output
  # out$random_effect_b <- list(postMeans = out_jags$mean$b,
  #                             postSd = out_jags$sd$b,
  #                             Rhat = out_jags$Rhat$b)
  # colnames(out$random_effect_b$postMeans) <-
  #   colnames(out$random_effect_b$postSd) <-
  #   colnames(out$random_effect_b$Rhat) <- colnames(Z1)
  #
  # # u random effect output
  # out$random_effect_u <- list(postMeans = out_jags$mean$u,
  #                             postSd = out_jags$sd$u,
  #                             Rhat = out_jags$Rhat$u)
  # colnames(out$random_effect_u$postMeans) <-
  #   colnames(out$random_effect_u$postSd) <-
  #   colnames(out$random_effect_u$Rhat) <- colnames(Z2)
  #
  # # median : Posterior median of parameters
  # out$median <- out_jags$q50
  # out$median$b <- NULL
  # out$median$u <- NULL
  #
  # # mean : Posterior mean of parameters
  # out$mean <- out_jags$mean
  # out$mean$b <- NULL
  # out$mean$u <- NULL
  #
  # # modes of parameters
  # out$modes <- lapply(out$sims.list, function(x) {
  #   m <- function(x) {
  #     d <- stats::density(x, bw = "nrd", adjust = 3, n = 1000)
  #     d$x[which.max(d$y)]
  #   }
  #   if (is.matrix(x))
  #     as.array(apply(x, 2, m))
  #   else{
  #     if(is.array(x))
  #       apply(x, c(2,3), m)
  #     else m(x)
  #   }
  # })
  #
  # # standard error of parameters
  # out$StErr <- lapply(out$sims.list, function(x) {
  #   f <- function(x) {
  #     acf.x <- drop(stats::acf(x, lag.max = 0.4 * length(x), plot = FALSE)$acf)[-1]
  #     acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
  #     ess <- length(x)/(1 + 2 * sum(acf.x))
  #     sqrt(stats::var(x)/ess)
  #   }
  #   if (is.matrix(x))
  #     as.array(apply(x, 2, f))
  #   else{
  #     if(is.array(x))
  #       apply(x, c(2,3), f)
  #     else f(x)
  #   }
  # })
  #
  # # standard deviation of parameters
  # out$StDev <- out_jags$sd
  # out$StDev$b <- NULL
  # out$StDev$u <- NULL
  #
  # # Rhat : Gelman & Rubin diagnostic
  # out$Rhat <- out_jags$Rhat
  # out$Rhat$b <- NULL
  # out$Rhat$u <- NULL
  #
  # # names
  # names(out$mean$beta) <-
  #   names(out$median$beta) <-
  #   names(out$modes$beta) <-
  #   names(out$StErr$beta) <-
  #   names(out$Rhat$beta) <-
  #   names(out$StDev$beta) <- colnames(X1)
  #
  # names(out$mean$xi) <-
  #   names(out$median$xi) <-
  #   names(out$modes$xi) <-
  #   names(out$StErr$xi) <-
  #   names(out$Rhat$xi) <-
  #   names(out$StDev$xi) <- colnames(X2)
  #
  # if(ncZ1==1){
  #   names(out$mean$covariance.b) <-
  #     names(out$median$covariance.b) <-
  #     names(out$modes$covariance.b) <-
  #     names(out$StErr$covariance.b) <-
  #     names(out$Rhat$covariance.b) <-
  #     names(out$StDev$covariance.b) <- colnames(Z1)
  # }else{
  #   colnames(out$mean$covariance.b) <-
  #     rownames(out$mean$covariance.b) <-
  #     colnames(out$median$covariance.b) <-
  #     rownames(out$median$covariance.b) <-
  #     colnames(out$modes$covariance.b) <-
  #     rownames(out$modes$covariance.b) <-
  #     colnames(out$Rhat$covariance.b) <-
  #     rownames(out$Rhat$covariance.b) <-
  #     colnames(out$StErr$covariance.b) <-
  #     rownames(out$StErr$covariance.b) <-
  #     colnames(out$StDev$covariance.b) <-
  #     rownames(out$StDev$covariance.b) <- colnames(Z1)
  # }
  #
  # if(ncZ2==1){
  #   names(out$mean$covariance.u) <-
  #     names(out$median$covariance.u) <-
  #     names(out$modes$covariance.u) <-
  #     names(out$StErr$covariance.u) <-
  #     names(out$Rhat$covariance.u) <-
  #     names(out$StDev$covariance.u) <- colnames(Z2)
  # }else{
  #   colnames(out$mean$covariance.u) <-
  #     rownames(out$mean$covariance.u) <-
  #     colnames(out$median$covariance.u) <-
  #     rownames(out$median$covariance.u) <-
  #     colnames(out$modes$covariance.u) <-
  #     rownames(out$modes$covariance.u) <-
  #     colnames(out$Rhat$covariance.u) <-
  #     rownames(out$Rhat$covariance.u) <-
  #     colnames(out$StErr$covariance.u) <-
  #     rownames(out$StErr$covariance.u) <-
  #     colnames(out$StDev$covariance.u) <-
  #     rownames(out$StDev$covariance.u) <- colnames(Z2)
  # }
  #
  # # credible intervals
  #
  # ## beta parameters
  # out$CIs$beta <- cbind(as.vector(t(out_jags$q2.5$beta)),
  #                       as.vector(t(out_jags$q97.5$beta)))
  # rownames(out$CIs$beta) <- colnames(X1)
  # colnames(out$CIs$beta) <- c("2.5%", "97.5%")
  #
  # ## xi parameters
  # out$CIs$xi <- cbind(as.vector(t(out_jags$q2.5$xi)),
  #                     as.vector(t(out_jags$q97.5$xi)))
  # rownames(out$CIs$xi) <- colnames(X2)
  # colnames(out$CIs$xi) <- c("2.5%", "97.5%")
  #
  # # only for diagonal elements of covariance matrix of random effects
  # if(ncZ1==1){
  #   out$CIs$covariance.b <- cbind(as.vector(t(out_jags$q2.5$covariance.b)),
  #                                 as.vector(t(out_jags$q97.5$covariance.b)))
  #   rownames(out$CIs$covariance.b) <- colnames(Z1)
  #   colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  # }else{
  #   out$CIs$covariance.b <- cbind(as.vector(diag(out_jags$q2.5$covariance.b)),
  #                                 as.vector(diag(out_jags$q97.5$covariance.b)))
  #   rownames(out$CIs$covariance.b) <- colnames(Z1)
  #   colnames(out$CIs$covariance.b) <- c("2.5%", "97.5%")
  # }
  # if(ncZ2==1){
  #   out$CIs$covariance.u <- cbind(as.vector(t(out_jags$q2.5$covariance.u)),
  #                                 as.vector(t(out_jags$q97.5$covariance.u)))
  #   rownames(out$CIs$covariance.u) <- colnames(Z2)
  #   colnames(out$CIs$covariance.u) <- c("2.5%", "97.5%")
  # }else{
  #   out$CIs$covariance.u <- cbind(as.vector(diag(out_jags$q2.5$covariance.u)),
  #                                 as.vector(diag(out_jags$q97.5$covariance.u)))
  #   rownames(out$CIs$covariance.u) <- colnames(Z2)
  #   colnames(out$CIs$covariance.u) <- c("2.5%", "97.5%")
  # }

  #---- End of the function defining the class and retruning the output
  class(out) <- c("Blsmm", "Bdr", "Bmm")
  out

}
