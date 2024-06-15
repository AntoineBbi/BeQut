#' @export
#'
summary.Blslqmm <- function (object, ...)
{
  if (!inherits(object, "Blslqmm"))
    stop("Use only with \"Blslqmm\" objects.\n")

  #---- Global details
  cat("#-- Statistical model: Location-scale linear quantile mixed model \n")
  cat(paste("     - Quantile order: ", object$control$tau, "\n"))
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
  cat(paste("     - Number of statistic units (e.g. subject): ", object$control$I, "\n"))
  cat("\n")

  #---- Parameter Estimations
  coefs <- object$mean
  CIs <- object$CIs
  Rhat <- object$Rhat
  # beta regression parameters
  beta_estim <- cbind(
    "Value" = coefs$beta,
    "2.5%" = CIs$beta[, 1],
    "97.5%" = CIs$beta[, 2],
    "Rhat" = Rhat$beta
  )
  sigma_estim <- cbind(
    "Value" = coefs$xi,
    "2.5%" = CIs$xi[, 1],
    "97.5%" = CIs$xi[, 2],
    "Rhat" = Rhat$xi
  )
  cat("#-- Estimation of location regression parameter(s) 'beta' and their CI bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat("#-- Estimation of scale regression parameter(s) 'xi': \n")
  prmatrix(sigma_estim, na.print = "")

  # Random effects for mixed regression model
  cat("\n")
  cat("#-- (Co)variance (matrix) of the location random-effect(s): \n")
  prmatrix(object$mean$covariance.b, na.print = "")

  cat("\n")
  cat("#-- (Co)variance (matrix) of the scale random-effect(s): \n")
  prmatrix(object$mean$covariance.u, na.print = "")

}
