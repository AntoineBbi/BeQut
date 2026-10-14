#' @export
#'
summary.Blqm <- function (object, ...)
{
  if (!inherits(object, "Blqm"))
    stop("Use only with \"Blqm\" objects.\n")

  #---- Global details
  cat("#-- Statistical model: Linear quantile model \n")
  cat(paste("     - Quantile order: ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
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
  cat("#-- Estimation of linear predictor parameters and their credible interval bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat(
    "#-- Estimation of 'sigma' parameter associated with asymmetric Laplace distribution: \n"
  )
  prmatrix(sigma_estim, na.print = "")

}
