#' @export
#'
predict.Blqmm <- function(object,
                          S = NULL, ...) {

  stopifnot(inherits(object, "Blqmm"))

  df     <- object$data
  id.var <- all.vars(object$control$formGroup)
  tau    <- object$control$tau
  post_mat   <- as.matrix(object$jags$out_jags$samples)
  frm    <- list(fixed = object$control$formFixed,
                 random = object$control$formRandom)

  # Get back objects to fit distribution parameters
  X <- object$jags$jags.data$X
  Z <- object$jags$jags.data$U
  offset <- object$jags$jags.data$offset
  I <- object$jags$jags.data$I
  y <- object$jags$jags.data$y
  N <- length(y)
  subj <- rep(seq_len(I), diff(offset))

  get_cols <- function(pattern) grep(pattern, colnames(post_mat))
  beta_cols  <- get_cols("^beta\\[")
  b_cols     <- get_cols("^b\\[")
  ncZ <- length(b_cols) / I

  # Get posterior mean of parameters and random effects and Compute distribution parameter for each observed Y
  beta_bar  <- colMeans(post_mat[, beta_cols,  drop = FALSE])
  b_bar <- matrix(colMeans(post_mat[, b_cols, drop = FALSE]), nrow = I, ncol = ncZ)
  mu   <- as.vector(X %*% beta_bar)  + rowSums(Z * b_bar[subj, , drop = FALSE])

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
  mu_draws <- matrix(NA_real_, nrow = S_eff, ncol = N)

  for (s in seq_len(S_eff)) {
    beta_s <- post_mat[s, beta_cols]
    b_s    <- matrix(post_mat[s, b_cols], nrow = I, ncol = ncZ)
    mu_draws[s, ] <- as.vector(X %*% beta_s) + rowSums(Z * b_s[subj, , drop = FALSE])
  }

  # output
  id <- df[[id.var]]
  pred <- data.frame(cbind(id = id,
                           y = y,
                           fitted.quantile = mu,
                           lower = apply(mu_draws, 2, quantile, probs = 0.025),
                           upper = apply(mu_draws, 2, quantile, probs = 0.975))
  )
  colnames(pred) <- c("id", "y", paste0("fitted.q", object$control$tau), "lower", "upper")
  pred
}
