#' @export
#'
summary.Blqmm <- function (object, ...)
{
  if (!inherits(object, "Blqmm"))
    stop("Use only with \"Blqmm\" objects.\n")

  #---- Global details
  cat("#-- Statistical model: Linear quantile mixed model \n")
  cat(paste("     - Quantile order (asymmetry parameter): ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
  cat(paste("     - Number of groups (e.g. individual): ", object$control$I, "\n"))
  cat("\n")

  # beta regression parameters
  beta_estim <- cbind(
    "Value" = as.vector(object$jags$out_jags$mean$beta),
    "2.5%" = as.vector(object$jags$out_jags$q2.5$beta),
    "97.5%" = as.vector(object$jags$out_jags$q97.5$beta),
    "Rhat" = as.vector(object$jags$out_jags$Rhat$beta)
  )
  rownames(beta_estim) <- colnames(object$jags$jags.data$X)
  sigma_estim <- cbind(
    "Value" = as.vector(object$jags$out_jags$mean$sigma),
    "2.5%" = as.vector(object$jags$out_jags$q2.5$sigma),
    "97.5%" = as.vector(object$jags$out_jags$q97.5$sigma),
    "Rhat" = as.vector(object$jags$out_jags$Rhat$sigma)
  )
  rownames(sigma_estim) <- "sigma"
  cat("#-- Estimation of longitudinal regression parameters and their credible interval bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat("#-- Estimation of scale parameter associated with asymmetric Laplace distribution: \n")
  prmatrix(sigma_estim, na.print = "")

  # Random effects for mixed regression model
  cat("\n")
  cat("#-- (Co)variance matrix of the random-effect(s): \n")
  if (object$control$RE_ind & is.matrix(object$jags$out_jags$mean$covariance.b)) {
    tmp <- diag(object$jags$out_jags$mean$covariance.b)
    colnames(tmp) <-
      rownames(tmp) <- colnames(object$jags$jags.data$U)
    prmatrix(tmp, na.print = "")
  }else{
    tmp <- object$jags$out_jags$mean$covariance.b
    colnames(tmp) <-
      rownames(tmp) <- colnames(object$jags$jags.data$U)
    prmatrix(tmp, na.print = "")
  }

}
