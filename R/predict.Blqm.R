#' @export
#'
predict.Blqm <- function(object,
                         newdata = NULL,
                         level = 0.95,
                         S = NULL, ...) {

  stopifnot(inherits(object, "Blqm"))

  if(!is.null(newdata) && !all(all.vars(object$control$formula)[-1] %in% colnames(newdata))){
     stop("'newdata' does not include all covariates used in the estimated model.\n")

  }
  if(!all.vars(object$control$formula)[1] %in% colnames(newdata)){
    newdata$y <- NA
    colnames(newdata)[ncol(newdata)] <- all.vars(object$control$formula)[1]
  }

  post_mat   <- as.matrix(object$jags$out_jags$samples)
  frm    <- object$control$formula

  # Get back objects to fit distribution parameters
  if(is.null(newdata)){
    X <- object$jags$jags.data$X
    I <- object$jags$jags.data$I
    y <- object$jags$jags.data$y
  }else{
    # data
    tmp <- newdata[c(all.vars(frm))]
    tmp <- as.data.frame(unique(tmp))
    y <- tmp[all.vars(frm)][, 1]
    I <- nrow(tmp)
    # design matrices
    frm2 <- frm
    frm2[[2]] <- NULL
    tmp <- as.data.frame(
      unique(newdata[all.vars(frm2)])
    )
    X <- stats::model.matrix(frm2, data = tmp)
  }

  get_cols <- function(pattern) grep(pattern, colnames(post_mat))
  beta_cols  <- get_cols("^beta\\[")

  # Get posterior mean of parameters and Compute distribution parameter for each observed Y
  beta_bar  <- colMeans(post_mat[, beta_cols,  drop = FALSE])
  mu   <- as.vector(X %*% beta_bar)

  ## draw iterations from posterior sample
  S_total <- nrow(post_mat)
  if (is.null(S)) {
    idx <- seq_len(S_total)
  } else {
    S   <- min(as.integer(S), S_total)
    idx <- round(seq(1, S_total, length.out = S))
  }
  post_mat <- post_mat[idx, , drop = FALSE]
  S_eff <- nrow(post_mat)

  ## conditionnal linear predictor
  mu_draws <- matrix(NA_real_, nrow = S_eff, ncol = I)

  for (s in seq_len(S_eff)) {
    beta_s <- post_mat[s, beta_cols]
    mu_draws[s, ] <- as.vector(X %*% beta_s)
  }

  # output
  alpha <- (1 - level) / 2
  pred <- data.frame(cbind(y = y,
                           fitted.quantile = mu,
                           lower = apply(mu_draws, 2, stats::quantile, probs = alpha),
                           upper = apply(mu_draws, 2, stats::quantile, probs = 1-alpha))
  )
  colnames(pred) <- c("y", paste0("fitted.q", object$control$tau), "lower", "upper")
  pred
}
