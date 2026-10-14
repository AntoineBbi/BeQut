#' @export
#'
ranef <- function(object, ...){

  if(!inherits(object, "Bmm")){
      stop("Use only with \"Bmm\" objects.\n")
  }else{
    # location component
    location <- list(postMeans = object$jags$out_jags$mean$b,
                     postSd = object$jags$out_jags$sd$b,
                     Rhat = object$jags$out_jags$Rhat$b)
    if(inherits(object, "Bdr")){
      out <- list(location = location)
      colnames(location$postMeans) <-
        colnames(location$postSd) <-
        colnames(location$Rhat) <- colnames(object$jags$jags.data$Z1)
      # scale component
      scale <- list(postMeans = object$jags$out_jags$mean$u,
                    postSd = object$jags$out_jags$sd$u,
                    Rhat = object$jags$out_jags$Rhat$u)
      colnames(scale$postMeans) <-
        colnames(scale$postSd) <-
        colnames(scale$Rhat) <- colnames(object$jags$jags.data$Z2)
      out$scale <- scale
      if(inherits(object, "Baldrm")){
        # skewness component
        skewness <- list(postMeans = object$jags$out_jags$mean$a,
                         postSd = object$jags$out_jags$sd$a,
                         Rhat = object$jags$out_jags$Rhat$a)
        colnames(skewness$postMeans) <-
          colnames(skewness$postSd) <-
          colnames(skewness$Rhat) <- colnames(object$jags$jags.data$Z3)
        out$skewness <- skewness
      }
    }else{
      colnames(location$postMeans) <-
        colnames(location$postSd) <-
        colnames(location$Rhat) <- colnames(object$jags$jags.data$U)
      out <- location
    }
  }
  out
}
