#' \code{lqmm} fits linear quantile mixed model
#'
#' Function using 'JAGS' software to estimate the linear quantile mixed model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formFixed formula for fixed part of longitudinal submodel with response variable
#' @param formRandom formula for random part of longitudinal submodel without response variable
#' @param formGroup formula specifying the cluster variable (e.g. = ~ subject)
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param RE_ind Boolean denoting if the random effects are assumed independent ; default is \code{FALSE}
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param object_lqm Object returns by 'lqm' function of BeQuT package to initialize parameter chains, default is \code{NULL}
#' @param precision precision (inverse of variance) of prior distribution, default is 0.1
#' @param parallel see \code{jagsUI::jags()} function
#'
#'
#' @return A \code{Blqmm} object is a list with the following elements:
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
#' @references Marco Geraci and Matteo Bottai (2014).
#' \emph{Linear quantile mixed models}.
#' Statistics and Computing, 24(3):461-479. doi: 10.1007/s11222-013-9381-9.
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
#' fit <- lqmm(formFixed = SBP ~ age.visit65,
#'             formRandom = ~ age.visit65,
#'             formGroup = ~ ID,
#'             data = threeC,
#'             tau = 0.75,
#'             n.iter = 1000,
#'             n.burnin = 500)
#'
#' #---- Get the posterior statistics
#' fit$summary_table
#'
#' #---- Get the summary of estimated model
#' summary(fit)
#'
#' #---- Visualize the trace for beta parameters
#' jagsUI::traceplot(fit$jags$out_jags, parameters = "beta")
#'
#' #---- Visualization of fitted quantile and its credible interval
#' pred <- plot(fit, time.var = "age.visit65")
#' }
#'
lqmm <- function(formFixed,
                  formRandom,
                  formGroup,
                  data,
                  tau,
                  RE_ind = FALSE,
                  n.chains = 3,
                  n.iter = 10000,
                  n.burnin = 5000,
                  n.thin = 1,
                  n.adapt = NULL,
                  object_lqm = NULL,
                  precision = 0.1,
                  parallel = FALSE){


  #-- data management
  data_long <- data[unique(c(all.vars(formGroup),all.vars(formFixed),all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- stats::model.frame(formFixed, data = data_long)
  X <- stats::model.matrix(formFixed, mfX)
  mfU <- stats::model.frame(formRandom, data = data_long)
  U <- stats::model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data_long)))
    data_long <- cbind(data_long, id = id)
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  # Initialization of prior parameter
  if(!is.null(object_lqm) && formFixed!=object_lqm$control$formula)
    stop("The 'formFixed' formula does not match the one from 'object_lqm'.")
  if(!is.null(object_lqm) && formFixed==object_lqm$control$formula){
    priorMean.beta <- object_lqm$mean$beta
    priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
  }else{
    priorMean.beta <- rep(0, ncol(X))
    priorTau.beta <- diag(rep(precision, length(priorMean.beta)))
  }
  # initialization of chains
  initial.values <- list(b = matrix(0, ncol = ncol(U), nrow = I),
                         beta = priorMean.beta,
                         sigma = ifelse(test = is.null(object_lqm),
                                        yes =  1,
                                        no = object_lqm$mean$sigma) #lqmm sigma = tmp_model$scale
  )

  # list of data jags
  jags.data <- list(y = y,
                    X = X,
                    U = U,
                    tau = tau,
                    ncX = ncol(X),
                    ncU = ncol(U),
                    I = I,
                    offset = offset,
                    priorMean.beta = priorMean.beta,
                    priorTau.beta = priorTau.beta,
                    priorA.sigma = precision,
                    priorB.sigma = precision
  )

  if(jags.data$ncU==1)
    RE_ind <- TRUE
  if(RE_ind){
    jags.data <- c(jags.data,
                   list(priorA.Sigma2 = precision,
                        priorB.Sigma2 = precision)
    )
    # initial.values$prec.Sigma2 <- 1/VarCorr(tmp_model) #lqmm
    initial.values$prec.Sigma2 <- rep(precision, ncol(U))
  }else{
    jags.data <- c(jags.data,
                   list(priorR.Sigma2 = diag(rep(precision, ncol(U))),
                        priorK.Sigma2 = ncol(U),
                        mu0 = rep(0, ncol(U))
                   )
    )
    # initial.values$prec.Sigma2 <- diag(1/VarCorr(tmp_model)) #lqmm
    initial.values$prec.Sigma2 <- diag(rep(precision, ncol(U)))
  }

  #---- write jags model in txt from R function
  if(RE_ind){
    jags_code <- "model{
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
  }#end of i loop
  # priors for parameters
  for(rr in 1:ncU){
    prec.Sigma2[rr] ~ dgamma(priorA.Sigma2, priorB.Sigma2)
    covariance.b[rr] <- 1/prec.Sigma2[rr]
  }
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
}"
  }else{
    jags_code <- "model{
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
  }#end of i loop
  # priors for parameters
  prec.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b <- inverse(prec.Sigma2[, ])
  beta[1:ncX] ~ dmnorm(priorMean.beta[], priorTau.beta[, ])
  sigma ~ dgamma(priorA.sigma, priorB.sigma)
}"
  }

  # regression from X
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * X[j, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], X[j, 1:ncX])", rplc, jags_code, fixed = TRUE)
  # regression on random effects b
  rplc <- paste(paste("b[i, ", 1:jags.data$ncU, "] * U[j, ", 1:jags.data$ncU, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(b[i, 1:ncU], U[j, 1:ncU])", rplc, jags_code, fixed = TRUE)

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
  parms_to_save <- c("beta", "sigma", "b", "covariance.b")

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
                      RE_ind = RE_ind,
                      I = I)

  # save jags information
  out$jags <- list(jags_code = jags_code,
                   jags.data = jags.data,
                   out_jags = out_jags)

  # summary output
  parms_to_save2 <- c("beta", "covariance.b")
  out$summary_table <- MCMCvis::MCMCsummary(out$jags$out_jags$samples,
                                            params = parms_to_save2)

    #---- End of the function defining the class and returning the output
  class(out) <- c("Blqmm", "Bmm")
  out

}
