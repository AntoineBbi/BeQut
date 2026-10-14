#' @export
#'
summary.Bdr <- function (object, ...)
{
  if (!inherits(object, "Bdr"))
    stop("Use only with \"Bdr\" objects.\n")

  #---- Global details
  if (inherits(object, "Blslqmm")){
    cat("#-- Statistical model: Location-scale linear quantile mixed model \n")
    cat(paste("     - Quantile order: ", object$control$tau, "\n"))
  }
  if (inherits(object, "Baldrm"))
    cat("#-- Statistical model: asymmetric Laplace distributional regresion mixed model \n")
  if (inherits(object, "Blsmm"))
    cat("#-- Statistical model: Location-scale mixed model  \n")
  cat(paste("     - Number of observations: ", nrow(object$data), "\n"))
  cat(paste("     - Number of statistic units (e.g. subject): ", object$control$I, "\n"))
  cat("\n")

  # location regression parameters
  beta_estim <- cbind(
    "Value" = as.vector(object$jags$out_jags$mean$beta),
    "2.5%" = as.vector(object$jags$out_jags$q2.5$beta),
    "97.5%" = as.vector(object$jags$out_jags$q97.5$beta),
    "Rhat" = as.vector(object$jags$out_jags$Rhat$beta)
  )
  rownames(beta_estim) <- colnames(object$jags$jags.data$X1)
  # scale regression parameters
  xi_estim <- cbind(
    "Value" = as.vector(object$jags$out_jags$mean$xi),
    "2.5%" = as.vector(object$jags$out_jags$q2.5$xi),
    "97.5%" = as.vector(object$jags$out_jags$q97.5$xi),
    "Rhat" = as.vector(object$jags$out_jags$Rhat$xi)
  )
  rownames(xi_estim) <- colnames(object$jags$jags.data$X2)
  cat("#-- Estimation of location regression parameter(s) 'beta' and their CI bounds: \n")
  prmatrix(beta_estim, na.print = "")
  cat("\n")
  cat("#-- Estimation of scale regression parameter(s) 'xi': \n")
  prmatrix(xi_estim, na.print = "")

  if (inherits(object, "Baldrm")){
    # skewness regression parameters
    alpha_estim <- cbind(
      "Value" = as.vector(object$jags$out_jags$mean$alpha),
      "2.5%" = as.vector(object$jags$out_jags$q2.5$alpha),
      "97.5%" = as.vector(object$jags$out_jags$q97.5$alpha),
      "Rhat" = as.vector(object$jags$out_jags$Rhat$alpha)
    )
    rownames(alpha_estim) <- colnames(object$jags$jags.data$X3)
    cat("\n")
    cat("#-- Estimation of skewness regression parameter(s) 'alpha': \n")
    prmatrix(alpha_estim, na.print = "")
  }

  # Random effects for mixed regression model
  cat("\n")
  cat("#-- (Co)variance (matrix) of the location random-effect(s): \n")
  tmp <- object$jags$out_jags$mean$covariance.b
  colnames(tmp) <-
    rownames(tmp) <- colnames(object$jags$jags.data$Z1)
  prmatrix(tmp, na.print = "")

  cat("\n")
  cat("#-- (Co)variance (matrix) of the scale random-effect(s): \n")
  tmp <- object$jags$out_jags$mean$covariance.u
  colnames(tmp) <-
    rownames(tmp) <- colnames(object$jags$jags.data$Z2)
  prmatrix(tmp, na.print = "")

  if (inherits(object, "Baldrm")){
    cat("\n")
    cat("#-- (Co)variance (matrix) of the skewness random-effect(s): \n")
    tmp <- object$jags$out_jags$mean$covariance.a
    colnames(tmp) <-
      rownames(tmp) <- colnames(object$jags$jags.data$Z3)
    prmatrix(tmp, na.print = "")
  }

}
