#' Plot posterior predictions from a Bayesian linear quantile regression
#'
#' @param object An object of class "Blqm".
#' @param var.x Name of the covariate to display on the x-axis.
#' @param n.grid Number of points used to construct the prediction grid.
#' @param S Number of posterior draws used for prediction.
#'   If NULL, all available iterations are used.
#' @param lwd Line width.
#' @param level Credibility level of the posterior interval.
#' @param ... Additional graphical arguments passed to geom_point().
#'
#' @return Invisibly returns a data frame containing posterior predictions.
#' @export

plot.Blqm <- function(object,
                      var.x,
                      n.grid = 100,
                      S = NULL,
                      lwd = 1,
                      level = 0.95,
                      ...) {

  stopifnot(inherits(object, "Blqm"))

  if(missing(var.x))
    stop("The argument 'var.x' must be specified.")
  if(length(var.x) != 1L || !is.character(var.x))
    stop("'var.x' must be a single variable name.")

  if(!is.numeric(n.grid) || length(n.grid) != 1L ||
      is.na(n.grid) || n.grid < 2 || n.grid != as.integer(n.grid))
    stop("'n.grid' must be an integer greater than or equal to 2.")

  if(!is.numeric(level) || length(level) != 1L ||
      is.na(level) || level <= 0 || level >= 1)
    stop("'level' must be between 0 and 1.")

  df   <- object$data
  post <- as.matrix(object$jags$out_jags$samples)
  frm  <- object$control$formula
  if(!var.x %in% names(df))
    stop(sprintf("'%s' is not a data variable.", var.x))

  # The covariate must belong to the right-hand side
  rhs <- stats::delete.response(stats::terms(frm))
  rhs_vars <- all.vars(rhs)
  if(!var.x %in% rhs_vars)
    stop("'var.x' must appear in the right-hand side of the formula.")


  # Select posterior iterations
  S_total <- nrow(post)
  if(is.null(S)){
    idx <- seq_len(S_total)
  }else{
    if (length(S) != 1L || !is.numeric(S) ||
        is.na(S) || !is.finite(S) ||
        S < 1 || S != as.integer(S))
      stop("'S' must be a positive integer or NULL.")
    S <- min(as.integer(S), S_total)
    idx <- unique(
      round(seq(1, S_total, length.out = S))
    )
  }
  post <- post[idx, , drop = FALSE]

  # Extract posterior draws of beta
  get_cols <- function(pattern) grep(pattern, colnames(post))
  beta_cols  <- get_cols("^beta\\[")

  # Get posterior draws of parameters
  beta_draws <- post[, beta_cols, drop = FALSE]

  # Construct the prediction grid
  y_name <- all.vars(frm[[2]])
  dsub <- df[, unique(c(y_name, rhs_vars)), drop = FALSE]
  x_min <- min(dsub[[var.x]])
  x_max <- max(dsub[[var.x]])
  x_grid <- seq(from = x_min,
                to = x_max,
                length.out = as.integer(n.grid)
                )

  # Fix all other covariates at representative values
  representative <- function(x) {
    # numeric variable
    if (is.numeric(x))
      return(stats::median(x, na.rm = TRUE))
    # factor variable
    if (is.factor(x)) {
      tab <- table(x)
      return(factor(
        names(tab)[which.max(tab)],
        levels = levels(x),
        ordered = is.ordered(x)
      ))
    }
    # logical variable
    if (is.logical(x))
      return(as.logical(names(which.max(table(x)))))

    x <- as.factor(x)
    tab <- table(x)
    factor(names(tab)[which.max(tab)], levels = levels(x))
  }


  other_vars <- setdiff(rhs_vars, var.x)
  ref_vals <- lapply(other_vars, function(v) representative(dsub[[v]]))
  names(ref_vals) <- other_vars

  newdata <- dsub[rep(1L, length(x_grid)), , drop = FALSE]
  for (v in other_vars) newdata[[v]] <- ref_vals[[v]]
  newdata[[var.x]] <- x_grid

  ref_txt <- if (length(other_vars) > 0) {
    paste0("Reference profile: ",
           paste(other_vars,
                 vapply(ref_vals, function(z) format(z, digits = 3), character(1)),
                 sep = " = ", collapse = ", "))
  } else NULL

  # Design matrix
  X <- stats::model.matrix(object$control$formula, data = newdata)


  # Posterior predictions
  # Rows = grid points; columns = posterior draws
  pred_draws <- X %*% t(beta_draws)
  pred_mean <- rowMeans(pred_draws)
  alpha <- (1 - level) / 2
  pred_low <- apply(pred_draws, 1, stats::quantile, probs = alpha, na.rm = TRUE)
  pred_high <- apply(pred_draws, 1, stats::quantile, probs = 1- alpha, na.rm = TRUE)

  pred <- data.frame(
    x = x_grid,
    fit = pred_mean,
    lower = pred_low,
    upper = pred_high
  )

  # plot

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dsub,
                        mapping = ggplot2::aes(x = .data[[var.x]], y = .data[[y_name]])) +
    ggplot2::geom_ribbon(data = pred,
                         mapping = ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         alpha = 0.2) +
    ggplot2::geom_line(data = pred,
                       mapping = ggplot2::aes(x = x, y = fit),
                       linewidth = lwd) +
    ggplot2::labs(x = var.x, y = y_name, caption = ref_txt) +
    ggplot2::theme_bw()

  print(p)

  invisible(pred)
}
