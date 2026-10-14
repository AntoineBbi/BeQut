#' \code{lqm} fits linear quantile regression model
#'
#' Function using 'JAGS' software to estimate the linear quantile regression model assuming asymmetric Laplace
#' distribution for residual error.
#'
#' @param formula formula for the quantile regression including response variable
#' @param data dataset of observed variables
#' @param tau the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non–duplicated values of \code{tau}, otherwise the execution is stopped.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 10000
#' @param n.burnin integer specifying how many of \code{n.iter} to discard as burn-in ; default is 5000
#' @param n.thin integer specifying the thinning of the chains; default is 1
#' @param n.adapt integer specifying the number of iterations to use for adaptation; default is \code{NULL}
#' @param parallel see \code{jagsUI::jags()} function
#'
#' @return A \code{Blqm} object which is a list with the following elements:
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
#' @export
#'
#' @examples
#'
#' \donttest{
#' #---- Use data
#' data(wave)
#'
#' #---- Fit regression model for the first quartile
#' lqm_025 <- BeQut::lqm(formula = h110d~vent_vit_moy + pression,
#'                       data = wave,
#'                       n.iter = 1000,
#'                       n.burnin = 500,
#'                       tau = 0.25)
#'
#' #---- Get the posterior statistics of parameters
#' lqm_025$summary_table
#'
#' #---- Plot the conditional first quartile
#' plot.Blqm(object = lqm_025, var.x = "vent_vit_moy")
#' plot.Blqm(object = lqm_025, var.x = "pression")
#' }
#'
lqm <- function(formula,
                data,
                tau = 0.5,
                n.chains = 3,
                n.iter = 10000,
                n.burnin = 5000,
                n.thin = 1,
                n.adapt = NULL,
                parallel = FALSE){

  # data
  tmp <- data[c(all.vars(formula))]
  tmp <- unique(tmp)
  I <- nrow(tmp)
  # design matrices
  y <- data[all.vars(formula)][,1]
  mfX <- stats::model.frame(formula, data = tmp)
  X <- stats::model.matrix(formula, mfX)
  ncX <- ncol(X)
  # use lm function to initiated values
  lm_tmp <- stats::lm(formula,
                      data = tmp)
  # prior beta parameters
  priorMean.beta <- as.numeric(lm_tmp$coefficients)
  # priorTau.beta <- diag(1/c(1,smcure_out$beta_var)/100)
  # list of data jags
  jags.data <- list(y = y,
                    X = X,
                    tau = tau,
                    ncX = ncX,
                    I = I
                    # ,
                    # priorMean.beta = priorMean.beta,
                    # priorTau.beta = priorTau.beta,
                    # priorA.shape = 1/100,
                    # priorB.shape = 1/100
                    )

  #---- write jags model in txt from R function
  jags_code <- "model{
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
}"
  rplc <- paste(paste("beta[", 1:jags.data$ncX, "] * X[i, ", 1:jags.data$ncX, "]", sep = ""), collapse = " + ")
  jags_code <- gsub("inprod(beta[1:ncX], X[i, 1:ncX])", rplc, jags_code, fixed = TRUE)

  # posterior samples to save
  parms_to_save <- c("beta", "sigma")

  # using jagsUI
  out_jags = jagsUI::jags(data = jags.data,
                          parameters.to.save = parms_to_save,
                          model.file = textConnection(jags_code),
                          n.chains = n.chains,
                          parallel = parallel,
                          n.adapt = n.adapt,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          DIC = F)

  #---- output
  out <- list(data = data)
  out$control <- list(formula = formula,
                      tau = tau,
                      call_function = "lqm",
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
  out$summary_table <- MCMCvis::MCMCsummary(out$jags$out_jags$samples,
                                            params = parms_to_save)

  #---- End of the function defining the class and returning the output
  class(out) <- c("Blqm")
  out

  }
