#' @export
#'
predict.Bdr <- function (object,
                         quantile = NULL,
                         expectation = FALSE,
                         variance = FALSE, ...){

  if (!inherits(object, "Bdr"))
    stop("Use only with \"Bdr\" objects.\n")

  # Get back objects to fit distribution parameters
  X1 <- object$jags$jags.data$X1
  Z1 <- object$jags$jags.data$Z1
  X2 <- object$jags$jags.data$X2
  Z2 <- object$jags$jags.data$Z2
  if(inherits(object, "Baldrm")){
    X3 <- object$jags$jags.data$X3
    Z3 <- object$jags$jags.data$Z3
  }
  offset <- object$jags$jags.data$offset
  I <- object$jags$jags.data$I
  y <- object$jags$jags.data$y
  N <- length(y)
  subj <- rep(seq_len(I), diff(offset))

  post_mat <- as.matrix(object$jags$out_jags$samples)
  get_cols <- function(pattern) grep(pattern, colnames(post_mat))
  beta_cols  <- get_cols("^beta\\[")
  xi_cols    <- get_cols("^xi\\[")
  b_cols     <- get_cols("^b\\[")
  u_cols     <- get_cols("^u\\[")
  ncZ1 <- length(b_cols) / I
  ncZ2 <- length(u_cols) / I
  if(inherits(object, "Baldrm")){
    alpha_cols <- get_cols("^alpha\\[")
    a_cols     <- get_cols("^a\\[")
    ncZ3 <- length(a_cols) / I
  }

  # Get posterior mean of parameters and random effects and Compute distribution parameter for each observed Y
  beta_bar  <- colMeans(post_mat[, beta_cols,  drop = FALSE])
  xi_bar    <- colMeans(post_mat[, xi_cols,    drop = FALSE])
  b_bar <- matrix(colMeans(post_mat[, b_cols, drop = FALSE]), nrow = I, ncol = ncZ1)
  u_bar <- matrix(colMeans(post_mat[, u_cols, drop = FALSE]), nrow = I, ncol = ncZ2)
  mu   <- as.vector(X1 %*% beta_bar)  + rowSums(Z1 * b_bar[subj, , drop = FALSE])
  sigma <- exp(as.vector(X2 %*% xi_bar) + rowSums(Z2 * u_bar[subj, , drop = FALSE]))
  if(inherits(object, "Baldrm")){
    alpha_bar <- colMeans(post_mat[, alpha_cols, drop = FALSE])
    a_bar <- matrix(colMeans(post_mat[, a_cols, drop = FALSE]), nrow = I, ncol = ncZ3)
    tau   <- plogis(as.vector(X3 %*% alpha_bar) + rowSums(Z3 * a_bar[subj, , drop = FALSE]))
  }
  if(inherits(object, "Blslqmm"))
    tau <- rep(object$control$tau, length(mu))


  if(is.numeric(quantile) && all(quantile > 0 & quantile < 1)){

    # when ALD is assumed
    if(inherits(object, "Bald")){
      # function computing quantile with order p from ALD(mu, sigma, tau)
      quantileALD <- function(mu, sigma, tau, p){
        Q_p <- ifelse(p<=tau,
                      (mu + sigma/(1-tau) * log(p/tau)),
                      (mu - sigma/tau * log((1-p)/(1-tau)))
        )
        return(Q_p)
      }

      fit_quantiles <- matrix(NA, nrow = length(mu), ncol = length(quantile))
      for(p in 1:length(quantile)){
        fit_quantiles[, p] <- quantileALD(mu, sigma, tau, p = quantile[p])
      }
    }

    # when Gaussian is assumed
    if(inherits(object, "Blsmm")){
      fit_quantiles <- matrix(NA, nrow = length(mu), ncol = length(quantile))
      for(p in 1:length(quantile)){
        fit_quantiles[, p] <- qnorm(p = quantile[p], mean = mu, sd = sigma)
      }
    }

    # add colnames to specify which quantile is predicted
    colnames(fit_quantiles) <- paste0("q_", quantile)
  }else{
    if(!is.null(quantile))
      stop("The ‘quantile’ argument must be either NULL or a vector of numeric values between 0 and 1 (exclusive),
           specifying the desired quantile orders.\n")
  }

  id <- object$data[all.vars(object$control$formGroup)][,1]
  out <- data.frame(cbind(id = id,
                          y = y,
                          mu = mu,
                          sigma = sigma)
                    )
  if(inherits(object, "Bald")){
    out$tau <- tau
    if(expectation)
      out$expectation <- out$mu + out$sigma*(1-2*out$tau)/(out$tau*(1-out$tau))
    if(variance)
      out$variance <- out$sigma^2*((out$tau^2+(1-out$tau)^2))/(out$tau^2*(1-out$tau)^2)
  }

  if(!is.null(quantile))
    out <- cbind(out, fit_quantiles)

  out

}
