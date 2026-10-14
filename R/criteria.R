#' \code{criteria} computes model criteria as DIC or WAIC from posterior samples.
#'
#'
#' @param object Object returns by \code{aldrm} function
#' @param name Character vector precising which criterion the user wants. By default, \code{name = c("DIC", "WAIC")} and both DIC and WAIC are computed. If \code{names} includes "calibration", the calibration criterion is returned given the quantile set and loss function chosen.
#' @param S Positive integer defining how much posterior draws is used to compute criteria. By default, \code{S = NULL} and the all posterior samples are used, but it can take time.
#' @param quantile Argument must be either NULL or a vector of numeric values between 0 and 1 (exclusive), specifying the desired quantile orders.
#' @param progress.bar By default is \code{TRUE}.
#'
#' @returns A list is a list with the following elements:
#'
#' @import dplyr
#' @importFrom loo waic
#'
#' @export
#'
#' @examples
#'
#' #' \donttest{
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
#' #---- Fit ALD regression model with subject-specific random effect defining all distribution parameters
#' fit <- aldrm(formFixed = SBP ~ age.visit65,
#'              formRandom = ~ age.visit65,
#'              formFixedScale = ~ age.visit65,
#'              formRandomScale = ~ age.visit65,
#'              formFixedSkewness = ~ age.visit65,
#'              formRandomSkewness = ~ age.visit65,
#'              formGroup = ~ ID,
#'              data = threeC,
#'              n.iter = 3000,
#'              n.burnin = 1000)
#'
#' #---- Fit criteria
#' crit_aldrm <- criteria(object = fit_aldrm,
#'                        name = c("DIC", "WAIC", "calibration"),
#'                        S = 500,
#'                        quantile = list(decile = seq(0.1,0.9, by=0.1),
#'                                        quartile = c(0.25,0.5,0.75)),
#'                        progress.bar = TRUE)
#'
#' }
#'
criteria <- function(object,
                     name = c("DIC", "WAIC", "calibration"),
                     S = NULL,
                     quantile = NULL,
                     progress.bar = TRUE) {

  if (!inherits(object, "Bdr"))
    stop("Use only with \"Bdr\" objects.\n")

  # initialization

  # out
  out <- list(control=list(object = object,
                           name = name,
                           quantile = quantile,
                           S = S))

  post_mat <- as.matrix(object$jags$out_jags$samples)
  S_total  <- nrow(post_mat)

  # subsample of iterations according S
  if (is.null(S)) {
    idx <- seq_len(S_total)
  } else {
    if (!is.numeric(S) || S < 0)
      stop("S must be an positive integer greater than 2.")
    S <- min(as.integer(S), S_total)
    idx <- round(seq(1, S_total, length.out = S))   # uniform sub-sample on the chain
  }
  post_mat <- post_mat[idx, , drop = FALSE]
  S_eff <- nrow(post_mat)

  # Get back objects for computing log-density
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

  if("DIC" %in% name || "WAIC" %in% name){
    # matrix/vector initialization
    ll_marg <- matrix(NA_real_, S_eff, N)
    if("DIC" %in% name)
      deviance_marg <- numeric(S_eff)

    pb <- if (progress.bar) txtProgressBar(min = 0, max = S_eff, style = 3, char = "-") else NULL

    # loop for each considered iteration to compute log-density
    for (s in seq_len(S_eff)) {

      # get current fe draws
      beta_s  <- post_mat[s, beta_cols]
      xi_s    <- post_mat[s, xi_cols]
      # get current re draws
      b_s <- matrix(post_mat[s, b_cols], nrow = I, ncol = ncZ1)
      u_s <- matrix(post_mat[s, u_cols], nrow = I, ncol = ncZ2)

      # compute current value of both location and scale parameters
      mu   <- as.vector(X1 %*% beta_s)  + rowSums(Z1 * b_s[subj, , drop = FALSE])
      sigma <- exp(as.vector(X2 %*% xi_s) + rowSums(Z2 * u_s[subj, , drop = FALSE]))

      # computation of Gaussian log-density
      if(inherits(object, "Blsmm")){
        resid2 <- (y - mu)^2
        sigma2 <- sigma^2
        ll_marg[s, ] <- -log(2*pi*sigma2)/2 - resid2 / (2*sigma2)
      }

      if(inherits(object, "Bald")){
        # compute current value of skewness parameter from current draws of fe and re

        if(inherits(object, "Baldrm")){
          alpha_s <- post_mat[s, alpha_cols]
          a_s <- matrix(post_mat[s, a_cols], nrow = I, ncol = ncZ3)
          tau   <- plogis(as.vector(X3 %*% alpha_s) + rowSums(Z3 * a_s[subj, , drop = FALSE]))
        }
        if(inherits(object, "Blslqmm"))
          tau <- object$control$tau
        # computation of log-density of ALD
        resid <- y - mu
        rho   <- resid * (tau - 1 * (resid < 0))
        ll_marg[s, ] <- log(tau) + log(1 - tau) - log(sigma) - rho / sigma
      }

      if("DIC" %in% name)
        deviance_marg[s] <- -2 * sum(ll_marg[s, ])

      if (progress.bar) setTxtProgressBar(pb, s)
    }
    if (progress.bar) close(pb)
  }


  if("DIC" %in% name){
    D_bar <- mean(deviance_marg, na.rm = TRUE)

    ## --- pD (Spiegelhalter): deviance with respect to theta_bar = posterior mean
    beta_bar  <- colMeans(post_mat[, beta_cols,  drop = FALSE])
    xi_bar    <- colMeans(post_mat[, xi_cols,    drop = FALSE])
    b_bar <- matrix(colMeans(post_mat[, b_cols, drop = FALSE]), nrow = I, ncol = ncZ1)
    u_bar <- matrix(colMeans(post_mat[, u_cols, drop = FALSE]), nrow = I, ncol = ncZ2)
    mu   <- as.vector(X1 %*% beta_bar)  + rowSums(Z1 * b_bar[subj, , drop = FALSE])
    sigma <- exp(as.vector(X2 %*% xi_bar) + rowSums(Z2 * u_bar[subj, , drop = FALSE]))

    # computation of gaussian log-density
    if(inherits(object, "Blsmm")){
      resid2 <- (y - mu)^2
      sigma2 <- sigma^2
      loglik_bar <- -log(2*pi*sigma2)/2 - resid2 / (2*sigma2)
    }

    if(inherits(object, "Bald")){
      if(inherits(object, "Baldrm")){
        alpha_bar <- colMeans(post_mat[, alpha_cols, drop = FALSE])
        a_bar <- matrix(colMeans(post_mat[, a_cols, drop = FALSE]), nrow = I, ncol = ncZ3)
        tau   <- plogis(as.vector(X3 %*% alpha_bar) + rowSums(Z3 * a_bar[subj, , drop = FALSE]))
      }
      if(inherits(object, "Blslqmm"))
        tau <- object$control$tau
      resid <- y - mu
      rho   <- resid * (tau - 1 * (resid < 0))
      loglik_bar <- log(tau) + log(1 - tau) - log(sigma) - rho / sigma
    }

    # computation of deviance with respect to posterior means
    D_theta_bar <- -2 * sum(loglik_bar)

    pD_spiegel <- D_bar - D_theta_bar
    pV_gelman  <- var(deviance_marg) / 2

    out$DIC <- list(DIC.spiegelhalter = D_bar + pD_spiegel,
                    DIC.gelman        = D_bar + pV_gelman,
                    pD.spiegelhalter  = pD_spiegel,
                    pD.gelman         = pV_gelman,
                    mean.deviance     = D_bar,
                    deviance.samples  = deviance_marg,
                    n.samples         = S_eff
    )
  }

  if("WAIC" %in% name){
    n_bad <- sum(!is.finite(ll_marg))
    if (n_bad > 0)
      warning(sprintf("loglik_marg: %d non finite values (NA/Inf) are detected.", n_bad))
    # compute waic of the model using waic function of loo package
    out$WAIC <- loo::waic(ll_marg)
  }


  # Calibration criteria
  if("calibration" %in% name && !is.null(quantile)){

    if(is.list(quantile)){

      # initialization
      n_set <- length(quantile)
      calibration <- vector("list", length = n_set)

      # computation for each set
      for(q in 1:n_set){
        if(is.numeric(quantile[[q]]) && all(quantile[[q]] > 0 & quantile[[q]] < 1)){

          vec <- sort(quantile[[q]])

          # prediction of quantile for each observations
          if(inherits(object, "Bdr"))
            pred <- predict.Bdr(object, quantile = vec)

          # compute the proportion of observation lower than the considered quantiles
          df_id <- pred %>%
            group_by(id) %>%
            summarise(
              n_measures = n(),
              across(
                starts_with("q_"),
                ~ sum(y < .x, na.rm = TRUE),
                .names = "n_y_inf_{.col}"
              ),
              .groups = "drop"
            )
          df_id <- as.data.frame(df_id)

          proportions <- df_id[, 3:(2+length(vec))] / df_id$n_measures
          fMAE <- function(v) mean(abs(v - vec))
          df_id$MAE <- apply(proportions, 1, fMAE)
          fMSE <- function(v) mean((v - vec)^2)
          df_id$MSE <- apply(proportions, 1, fMSE)

          MMAE <- mean(df_id$MAE)
          MMSE <- mean(df_id$MSE)

          out$calibration[[q]] <- list(df_id = df_id,
                                       MMAE = MMAE,
                                       MMSE = MMSE,
                                       quantile = vec)

        }else{
          if(!is.null(quantile))
            stop("The ‘quantile’ argument (or elements of 'quantile' if 'quantile' is a list of quantile set) must be either NULL or a vector of numeric values between 0 and 1 (exclusive),
           specifying the desired quantile orders.\n")
        }
      }

      if(is.null(names(quantile))){
        names(out$calibration) <- paste0("Set", 1:n_set)
      }else{
        names(out$calibration) <- names(quantile)
      }

    }else{

      if(is.numeric(quantile) && all(quantile > 0 & quantile < 1)){

        vec <- sort(quantile)

        # prediction of quantile for each observations
        if(inherits(object, "Baldrm"))
          pred <- predict.Bdr(object, quantile = vec)

        # compute the proportion of observation lower than the considered quantiles
        df_id <- pred %>%
          group_by(id) %>%
          summarise(
            n_measures = n(),
            across(
              starts_with("q_"),
              ~ sum(y < .x, na.rm = TRUE),
              .names = "n_y_inf_{.col}"
            ),
            .groups = "drop"
          )
        df_id <- as.data.frame(df_id)

        proportions <- df_id[, 3:(2+length(vec))] / df_id$n_measures
        fMAE <- function(v) mean(abs(v - vec))
        df_id$MAE <- apply(proportions, 1, fMAE)
        fMSE <- function(v) mean((v - vec)^2)
        df_id$MSE <- apply(proportions, 1, fMSE)

        MMAE <- mean(df_id$MAE)
        MMSE <- mean(df_id$MSE)

        out$calibration <- list(df_id = df_id,
                                MMAE = MMAE,
                                MMSE = MMSE,
                                quantile = vec)

      }else{
        if(!is.null(quantile))
          stop("The ‘quantile’ argument (or elements of 'quantile' if 'quantile' is a list of quantile set) must be either NULL or a vector of numeric values between 0 and 1 (exclusive),
           specifying the desired quantile orders.\n")
      }

    }

  }

  # return the output in list
  out
}
