#' @export
#'
plot.Blqmm <- function(object,
                       individuals = NULL,
                       time.var,
                       n.grid = 100,
                       S = NULL,
                       level = 0.95,
                       show.observed = TRUE,
                       lwd = 1, ...) {

  stopifnot(inherits(object, "Blqmm"))
  if (missing(time.var))
    stop("Time variable must be specified (argument 'time.var').")

  df     <- object$data
  id.var <- all.vars(object$control$formGroup)
  tau    <- object$control$tau
  post   <- as.matrix(object$jags$out_jags$samples)
  frm    <- list(fixed = object$control$formFixed,
                 random = object$control$formRandom)

  if (!time.var %in% names(df))
    stop(sprintf("'%s' is not a data variable.", time.var))

  id_levels <- (unique(df[[id.var]]))
  I <- length(id_levels)

  # id <- as.integer(df[[id.var]])

  ## by default, draw randomly 9 individuals
  if (is.null(individuals)) {
    individuals <- sample(id_levels, 9)
  } else if (!all(individuals %in% id_levels)) {
    stop("Some id in 'individuals' are not observed.")
  }

  ## draw iterations from posterior sample
  S_total <- nrow(post)
  if (is.null(S)) {
    idx <- seq_len(S_total)
  } else {
    S   <- min(as.integer(S), S_total)
    idx <- round(seq(1, S_total, length.out = S))
  }
  post <- post[idx, , drop = FALSE]
  S_eff <- nrow(post)

  ncZ <- ncol(model.matrix(frm$random, df))
  I <- length(unique(df[[id.var]]))
  beta_cols <- grep("^beta\\[", colnames(post))
  b_cols    <- grep("^b\\[",    colnames(post))

  # # posterior means for both fixed and random parameters
  # get_mean <- function(pattern) colMeans(post[, grep(pattern, colnames(post)), drop = FALSE])
  # beta_bar  <- get_mean("^beta\\[")
  # b_bar <- matrix(get_mean("^b\\["), nrow = I, ncol = ncZ)

  traj_list <- lapply(seq_along(individuals), function(k) {
    id_i  <- individuals[k]
    row_i <- match(id_i, id_levels)
    dsub  <- df[df[[id.var]] == id_i, , drop = FALSE]

    t_min <- min(dsub[[time.var]])
    t_max <- max(dsub[[time.var]])
    t_grid <- seq(t_min, t_max, length.out = n.grid)

    newdata <- dsub[rep(1, n.grid), , drop = FALSE]
    newdata[[time.var]] <- t_grid
    row.names(newdata) <- NULL

    X_t <- model.matrix(frm$fixed,  newdata)
    Z_t <- model.matrix(frm$random, newdata)
    # mu_t_marg    <- as.vector(X_t %*% beta_bar)
    # mu_t    <- as.vector(X_t %*% beta_bar)  + as.vector(Z_t %*% b_bar[row_i, ])

    ## marginal and conditionnal linear predictor
    mu_marg <- matrix(NA_real_, S_eff, n.grid)   # X'beta (population)
    mu_cond <- matrix(NA_real_, S_eff, n.grid)   # X'beta + Z'b_i (id)

    for (s in seq_len(S_eff)) {
      beta_s <- post[s, beta_cols]
      b_s    <- matrix(post[s, b_cols], nrow = I, ncol = ncZ)[row_i, ]

      mu_marg[s, ] <- as.vector(X_t %*% beta_s)
      mu_cond[s, ] <- mu_marg[s, ] + as.vector(Z_t %*% b_s)
    }


    alpha <- (1 - level) / 2

    ## summary
    summarise <- function(mat, label) {
      data.frame(
        id = id_i,
        time = t_grid,
        type = label,
        # fitted.y  = ifelse(label=="marginal", mu_t_marg, mu_t),
        fitted.y  = colMeans(mat),
        lower = apply(mat, 2, quantile, probs = alpha),
        upper = apply(mat, 2, quantile, probs = (1-alpha))
      )
    }

    rbind(summarise(mu_marg, "marginal"), summarise(mu_cond, "conditional"))
  })

  pred <- do.call(rbind, traj_list)
  attr(pred, "tau") <- tau
  class(pred) <- c("lqmm.predict", "data.frame")

  tau <- attr(pred, "tau")

  p <- ggplot2::ggplot(pred, ggplot2::aes(x = time,
                                          y = fitted.y,
                                          color = type,
                                          fill = type)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = lwd) +
    ggplot2::facet_wrap(~ id) +
    ggplot2::labs(x = "Time",
                  y = sprintf("Quantile with order %.2f (IC 95%%)", tau),
                  color = "Quantile", fill = "Quantile") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

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
  invisible(pred)
}
