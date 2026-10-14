#' @export
#'
plot.Bdr<- function(object,
                    quantiles = c(0.25, 0.5, 0.75),
                    individuals,
                    time.var,
                    n.grid = 100,
                    level = 0.95,
                    show.observed = TRUE,
                    show.mode = TRUE,
                    S = NULL,
                    lwd = 1, ...) {

  # object = fit_aldrm
  # quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)
  # show.observed = T
  # individuals = individuals
  # time.var = "time"
  # n.grid = 100
  # show.observed = TRUE
  # show.mode = T
  # S = NULL
  # lwd = 1


  stopifnot(inherits(object, "Bdr"))
  if (missing(individuals))
    stop("Individual must be specified (argument 'individuals').")
  if (missing(time.var))
    stop("Time variable must be specified (argument 'time.var').")

  df     <- object$data
  id.var <- all.vars(object$control$formGroup)
  post   <- as.matrix(object$jags$out_jags$samples)

  if (!time.var %in% names(df))
    stop(sprintf("'%s' is not include in data used to fit the model.", time.var))
  if (!all(individuals %in% unique(df[[id.var]])))
    stop("Some id in 'individuals' are missing.")

  ## --- Sous-echantillonnage des tirages (cf. compute_waic_marg) ------------
  if (!is.null(S)) {
    S_total <- nrow(post)
    S   <- min(as.integer(S), S_total)
    idx <- round(seq(1, S_total, length.out = S))
    post_mat <- post[idx, , drop = FALSE]
    S_eff <- nrow(post_mat)
    beta_cols  <- grep("^beta\\[",  colnames(post_mat))
    xi_cols    <- grep("^xi\\[",    colnames(post_mat))
    alpha_cols <- grep("^alpha\\[", colnames(post_mat))
    b_cols     <- grep("^b\\[",     colnames(post_mat))
    u_cols     <- grep("^u\\[",     colnames(post_mat))
    a_cols     <- grep("^a\\[",     colnames(post_mat))
  }

  # posterior means for both fixed and random parameters
  get_mean <- function(pattern) colMeans(post[, grep(pattern, colnames(post)), drop = FALSE])
  beta_bar  <- get_mean("^beta\\[")
  xi_bar    <- get_mean("^xi\\[")

  ncZ1 <- ncol(model.matrix(object$control$formRandom,    df))
  ncZ2 <- ncol(model.matrix(object$control$formRandomScale, df))
  I <- length(unique(df[[id.var]]))

  b_bar <- matrix(get_mean("^b\\["), nrow = I, ncol = ncZ1)
  u_bar <- matrix(get_mean("^u\\["), nrow = I, ncol = ncZ2)


  # aldrm
  if(inherits(object, "Baldrm")){
    alpha_bar <- get_mean("^alpha\\[")
    ncZ3 <- ncol(model.matrix(object$control$formRandomSkewness,   df))
    a_bar <- matrix(get_mean("^a\\["), nrow = I, ncol = ncZ3)
  }

  id_levels <- sort(unique(df[[id.var]]))


  if(inherits(object, "Bald"))
    # Quantile of ALD
    qAL <- function(p, mu, sigma, tau) {
      ifelse(p <= tau,
             mu + sigma / (1 - tau) * log(p / tau),
             mu - sigma / tau       * log((1 - p) / (1 - tau)))
    }

  # Build trajectories by individuals
  traj_list <- lapply(individuals, function(id_i) {

    row_i <- match(id_i, id_levels)
    dsub  <- df[df[[id.var]] == id_i, , drop = FALSE]

    t_min <- min(dsub[[time.var]])
    t_max <- max(dsub[[time.var]])
    t_grid <- seq(t_min, t_max, length.out = n.grid)

    newdata <- dsub[rep(1, n.grid), , drop = FALSE]
    newdata[[time.var]] <- t_grid
    row.names(newdata) <- NULL

    # matrices design from saved formula
    X1_t <- model.matrix(object$control$formFixed, newdata)
    X2_t <- model.matrix(object$control$formFixedScale, newdata)

    Z1_t <- model.matrix(object$control$formRandom, newdata)
    Z2_t <- model.matrix(object$control$formRandomScale, newdata)

    mu_t    <- as.vector(X1_t %*% beta_bar)  + as.vector(Z1_t %*% b_bar[row_i, ])
    sigma_t <- exp(as.vector(X2_t %*% xi_bar) + as.vector(Z2_t %*% u_bar[row_i, ]))

    if(inherits(object, "Baldrm")){
      X3_t <- model.matrix(object$control$formFixedSkewness, newdata)
      Z3_t <- model.matrix(object$control$formRandomSkewness, newdata)
      tau_t   <- plogis(as.vector(X3_t %*% alpha_bar) + as.vector(Z3_t %*% a_bar[row_i, ]))
    }
    if(inherits(object, "Blslqmm"))
      tau_t <- rep(object$control$tau, length(mu_t))

    if(inherits(object, "Bald")){
      q_mat <- sapply(quantiles, function(p) qAL(p, mu_t, sigma_t, tau_t))
      colnames(q_mat) <- paste0("q", quantiles)
    }
    if(inherits(object, "Blsmm")){
      q_mat <- sapply(quantiles, function(p) qnorm(p = p, mean = mu_t, sd = sigma_t))
      colnames(q_mat) <- paste0("q", quantiles)
    }

    if(is.null(S)){

      data.frame(id = id_i, time = t_grid, mode = mu_t, q_mat,
                 check.names = FALSE)

    }else{

      ## --- Tableaux (S_eff x n.grid) : mode + chaque quantile, par tirage ----
      mode_mat <- matrix(NA_real_, S_eff, n.grid)
      q_arr    <- array(NA_real_, dim = c(S_eff, n.grid, length(quantiles)))

      for (s in seq_len(S_eff)) {
        beta_s  <- post_mat[s, beta_cols]
        xi_s    <- post_mat[s, xi_cols]

        b_s <- matrix(post_mat[s, b_cols], nrow = I, ncol = ncZ1)[row_i, ]
        u_s <- matrix(post_mat[s, u_cols], nrow = I, ncol = ncZ2)[row_i, ]

        mu_s    <- as.vector(X1_t %*% beta_s)  + as.vector(Z1_t %*% b_s)
        sigma_s <- exp(as.vector(X2_t %*% xi_s) + as.vector(Z2_t %*% u_s))

        if(inherits(object, "Baldrm")){
          alpha_s <- post_mat[s, alpha_cols]
          a_s <- matrix(post_mat[s, a_cols], nrow = I, ncol = ncZ3)[row_i, ]
          tau_s   <- plogis(as.vector(X3_t %*% alpha_s) + as.vector(Z3_t %*% a_s))
        }
        if(inherits(object, "Blslqmm"))
          tau_s <- rep(object$control$tau, length(mu_t))

        mode_mat[s, ] <- mu_s

        if(inherits(object, "Bald")){
          for (qi in seq_along(quantiles))
            q_arr[s, , qi] <- qAL(quantiles[qi], mu_s, sigma_s, tau_s)
        }
        if(inherits(object, "Blsmm")){
          for (qi in seq_along(quantiles))
            q_arr[s, , qi] <- qnorm(p = quantiles[qi], mean = mu_s, sd = sigma_s)
        }

      }

      ## --- Resume : moyenne posterior + IC a 95% (quantiles empiriques) -----

      alpha <- (1 - level) / 2

      mode_summary <- data.frame(
        id = id_i,
        time = t_grid,
        quantile = NA_character_,
        mean  = mu_t,
        lower = apply(mode_mat, 2, quantile, probs = alpha),
        upper = apply(mode_mat, 2, quantile, probs = 1-alpha),
        is_mode = TRUE
      )

      q_summary <- do.call(rbind, lapply(seq_along(quantiles), function(qi) {
        data.frame(
          id = id_i,
          time = t_grid,
          quantile = paste0("q", quantiles[qi]),
          mean  = q_mat[ ,qi],
          lower = apply(q_arr[, , qi], 2, quantile, probs = alpha),
          upper = apply(q_arr[, , qi], 2, quantile, probs = 1-alpha),
          is_mode = FALSE
        )
      }))

      rbind(mode_summary, q_summary)
    }


  })

  traj_df <- do.call(rbind, traj_list)

  ## Long format for ggplot2
  if(is.null(S)){

    q_cols <- paste0("q", quantiles)
    traj_long <- reshape(traj_df,
                         direction = "long",
                         varying = q_cols,
                         v.names = "value",
                         timevar = "quantile",
                         times = quantiles,
                         idvar = c("id", "time"))

    p <- ggplot2::ggplot(traj_long, ggplot2::aes(x = time, y = value,
                                                 color = factor(quantile))) +
      ggplot2::geom_line(lwd = lwd) +
      viridis::scale_color_viridis(discrete = TRUE) +
      ggplot2::geom_line(data = unique(traj_df[, c("id", "time", "mode")]),
                         ggplot2::aes(x = time,
                                      y = mode),
                         inherit.aes = FALSE,
                         color = "black",
                         linetype = "dashed", lwd = lwd) +
      ggplot2::facet_wrap(~ id) +
      ggplot2::labs(x = "Time",
                    y = "Response Y",
                    color = "Quantile'\n order") +
      ggplot2::theme_bw()

    if (show.observed) {
      obs <- object$data[object$data[[id.var]] %in% individuals, ]
      obs$id <- obs[[id.var]]
      y.var <- all.vars(object$control$formFixed)[1]
      p <- p + ggplot2::geom_point(data = obs,
                                   ggplot2::aes(x = .data[[time.var]],
                                                y = .data[[y.var]]),
                                   inherit.aes = FALSE, alpha = 0.3, size = 0.8)

    }

    print(p)
    invisible(traj_df)
  }else{
    q_data <- traj_df[!traj_df$is_mode, ]
    mode_data <- traj_df[traj_df$is_mode, ]

    p <- ggplot2::ggplot(q_data, ggplot2::aes(x = time,
                                              y = mean,
                                              color = quantile,
                                              fill = quantile)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = 0.15, color = NA) +
      ggplot2::geom_line(linewidth = lwd) +
      viridis::scale_color_viridis(discrete = TRUE) +
      viridis::scale_fill_viridis(discrete = TRUE) +
      # ggplot2::facet_grid(id ~ quantile) +
      ggplot2::facet_wrap(~ id) +
      # ggplot2::labs(x = "Time",
      #               y = "Response Y",
      #               color = "Quantile'\n order") +
      ggplot2::labs(x = "Time",
                    y = "Fitted value with CI at 95%",
                    color = "Quantile",
                    fill = "Quantile") +
      # ggplot2::theme(legend.position = "none") +
      ggplot2::theme_bw()

    if (show.mode) {
      ## le mode est identique pour chaque facette "quantile" d'un meme id ;
      ## on le duplique sur la grille de facettes pour l'afficher partout
      mode_rep <- do.call(rbind, lapply(paste0("q", quantiles), function(q) {
        d <- mode_data
        d$quantile <- q
        d
      }))
      p <- p + ggplot2::geom_line(data = mode_rep,
                                  ggplot2::aes(x = time, y = mean),
                                  inherit.aes = FALSE, color = "black",
                                  linetype = "dashed",
                                  linewidth = lwd)
    }

    if (show.observed) {

      obs <- object$data[object$data[[id.var]] %in% individuals, ]
      obs$id <- obs[[id.var]]
      y.var <- all.vars(object$control$formFixed)[1]
      p <- p + ggplot2::geom_point(data = obs,
                                   ggplot2::aes(x = .data[[time.var]],
                                                y = .data[[y.var]]),
                                   inherit.aes = FALSE, alpha = 0.3, size = 0.8)

    }

    print(p)
    invisible(traj_df)
  }

}
