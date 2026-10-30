# Extreme QR competitor: arm-wise IPW quantile regression, with extremal subsampling intervals
# (Chernozhukov and Fernandez-Val, 2011) at extreme levels and bootstrap intervals otherwise.
extremeqr_build_arm_formula <- function(y_name, x_names) {
  rhs <- if (length(x_names)) paste(x_names, collapse = " + ") else "1"
  as.formula(paste0(y_name, " ~ ", rhs))
}

extremalqr_should_use_extremal <- function(tau, n, d, threshold = 20, tail_cutoff = 0.2) {
  tail_prob <- min(tau, 1 - tau)
  order_adj <- tail_prob * n / max(d, 1)
  is.finite(order_adj) && tail_prob < tail_cutoff && order_adj <= threshold
}

extremalqr_subsample_size <- function(n, d) {
  min(n - 1L, max(d + 10L, floor(0.2 * n)))
}

extremalqr_weighted_mean <- function(x, weights = NULL) {
  if (is.null(weights)) {
    return(colMeans(x))
  }
  w <- as.numeric(weights)
  if (length(w) != nrow(x) || !all(is.finite(w)) || sum(w) <= 0) {
    return(rep(NA_real_, ncol(x)))
  }
  drop(crossprod(w, x) / sum(w))
}

extremalqr_compute_ps_weights <- function(d, x, trim = 1e-4, cap = NULL, stabilized = TRUE) {
  x_df <- as.data.frame(x)
  if (is.null(names(x_df))) {
    names(x_df) <- paste0("X", seq_len(ncol(x_df)))
  }

  ps_fit <- tryCatch(
    stats::glm(D ~ ., data = data.frame(D = d, x_df, check.names = FALSE), family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(ps_fit)) {
    return(rep(1, length(d)))
  }

  ehat <- as.numeric(stats::predict(ps_fit, type = "response"))
  ehat <- pmin(pmax(ehat, trim), 1 - trim)
  p1 <- mean(d == 1)

  weights <- if (stabilized) {
    ifelse(d == 1, p1 / ehat, (1 - p1) / (1 - ehat))
  } else {
    ifelse(d == 1, 1 / ehat, 1 / (1 - ehat))
  }

  if (!is.null(cap)) {
    weights <- pmin(weights, cap)
  }

  weights / mean(weights)
}

extremalqr_fit_beta <- function(y, x, tau, weights = NULL) {
  tau <- min(max(as.numeric(tau), 1e-4), 0.999)
  fit <- tryCatch(
    if (is.null(weights)) {
      quantreg::rq.fit.br(x = x, y = y, tau = tau)
    } else {
      quantreg::rq.wfit(x = x, y = y, tau = tau, weights = weights, method = "br")
    },
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$coef)) {
    return(NULL)
  }
  beta <- as.numeric(fit$coef)
  if (length(beta) != ncol(x) || any(!is.finite(beta))) {
    return(NULL)
  }
  beta
}

extremeqr_prepare_common_model <- function(x) {
  x_df <- as.data.frame(x)
  if (is.null(names(x_df))) {
    names(x_df) <- paste0("X", seq_len(ncol(x_df)))
  }
  train_df <- data.frame(Y = 0, x_df, check.names = FALSE)
  fit_formula <- extremeqr_build_arm_formula("Y", names(x_df))
  terms_x <- stats::delete.response(stats::terms(fit_formula, data = train_df))
  list(x_df = x_df, terms_x = terms_x)
}

extremeqr_model_matrix <- function(terms_x, x_df) {
  stats::model.matrix(terms_x, data = as.data.frame(x_df, check.names = FALSE))
}

extremeqr_align_new_x <- function(x_df, new_x) {
  new_x_df <- as.data.frame(new_x, check.names = FALSE)
  if (!identical(colnames(x_df), colnames(new_x_df))) {
    colnames(new_x_df) <- colnames(x_df)
  }
  new_x_df
}

extremeqr_split_weights <- function(weights, d_vec) {
  if (is.null(weights)) {
    return(list(treat = NULL, control = NULL))
  }
  list(
    treat = weights[d_vec == 1],
    control = weights[d_vec == 0]
  )
}

extremeqr_arm_boot_interval <- function(y, d, x, new_x, tau, alpha = 0.1,
                                        n_boot = 200, seed = NULL,
                                        weights = NULL, reestimate_ps = FALSE,
                                        ps_trim = 1e-4, ps_cap = NULL) {
  prep <- extremeqr_prepare_common_model(x)
  x_df <- prep$x_df
  new_x_df <- extremeqr_align_new_x(x_df, new_x)
  phi_mat <- extremeqr_model_matrix(prep$terms_x, new_x_df)

  idx_treat <- d == 1
  idx_control <- d == 0
  if (!any(idx_treat) || !any(idx_control)) {
    return(NULL)
  }

  x_treat <- extremeqr_model_matrix(prep$terms_x, x_df[idx_treat, , drop = FALSE])
  x_control <- extremeqr_model_matrix(prep$terms_x, x_df[idx_control, , drop = FALSE])
  arm_weights <- extremeqr_split_weights(weights, d)

  beta_treat <- extremalqr_fit_beta(y[idx_treat], x_treat, tau, weights = arm_weights$treat)
  beta_control <- extremalqr_fit_beta(y[idx_control], x_control, tau, weights = arm_weights$control)
  if (is.null(beta_treat) || is.null(beta_control)) {
    return(NULL)
  }

  theta_hat <- as.numeric(phi_mat %*% beta_treat) - as.numeric(phi_mat %*% beta_control)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- nrow(x_df)
  boot_draws <- matrix(NA_real_, nrow = nrow(phi_mat), ncol = n_boot)

  for (b_idx in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    y_boot <- y[idx]
    d_boot <- d[idx]
    x_boot_df <- x_df[idx, , drop = FALSE]

    if (length(unique(d_boot)) < 2L) {
      next
    }

    weights_boot <- NULL
    if (reestimate_ps) {
      weights_boot <- extremalqr_compute_ps_weights(
        d = d_boot,
        x = x_boot_df,
        trim = ps_trim,
        cap = ps_cap
      )
    } else if (!is.null(weights)) {
      weights_boot <- weights[idx]
    }
    arm_weights_boot <- extremeqr_split_weights(weights_boot, d_boot)

    idx_treat_boot <- d_boot == 1
    idx_control_boot <- d_boot == 0
    x_treat_boot <- tryCatch(
      extremeqr_model_matrix(prep$terms_x, x_boot_df[idx_treat_boot, , drop = FALSE]),
      error = function(e) NULL
    )
    x_control_boot <- tryCatch(
      extremeqr_model_matrix(prep$terms_x, x_boot_df[idx_control_boot, , drop = FALSE]),
      error = function(e) NULL
    )
    if (is.null(x_treat_boot) || is.null(x_control_boot)) {
      next
    }

    beta_treat_boot <- extremalqr_fit_beta(
      y_boot[idx_treat_boot], x_treat_boot, tau, weights = arm_weights_boot$treat
    )
    beta_control_boot <- extremalqr_fit_beta(
      y_boot[idx_control_boot], x_control_boot, tau, weights = arm_weights_boot$control
    )
    if (is.null(beta_treat_boot) || is.null(beta_control_boot)) {
      next
    }

    boot_draws[, b_idx] <-
      as.numeric(phi_mat %*% beta_treat_boot) -
      as.numeric(phi_mat %*% beta_control_boot)
  }

  valid_cols <- colSums(is.finite(boot_draws)) == nrow(boot_draws)
  boot_draws <- boot_draws[, valid_cols, drop = FALSE]
  if (ncol(boot_draws) < max(25, ceiling(n_boot / 5))) {
    return(NULL)
  }

  data.frame(
    estimate = theta_hat,
    sd = apply(boot_draws, 1, stats::sd, na.rm = TRUE),
    lower = apply(boot_draws, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE),
    upper = apply(boot_draws, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  )
}

extremeqr_arm_sn_interval <- function(y, d, x, new_x, tau, alpha = 0.1,
                                      spacing_p = 5, b = NULL,
                                      n_subsamples = 200, seed = NULL,
                                      weights = NULL, reestimate_ps = FALSE,
                                      ps_trim = 1e-4, ps_cap = NULL) {
  prep <- extremeqr_prepare_common_model(x)
  x_df <- prep$x_df
  new_x_df <- extremeqr_align_new_x(x_df, new_x)
  phi_mat <- extremeqr_model_matrix(prep$terms_x, new_x_df)
  x_mat_full <- extremeqr_model_matrix(prep$terms_x, x_df)

  n <- nrow(x_mat_full)
  d_coef <- ncol(x_mat_full)
  if (is.null(b)) {
    b <- extremalqr_subsample_size(n, d_coef)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  mtau <- tau + (d_coef + spacing_p) / n
  if (tau <= 0 || tau >= 0.2 || mtau >= 1) {
    return(NULL)
  }

  idx_treat <- d == 1
  idx_control <- d == 0
  if (!any(idx_treat) || !any(idx_control)) {
    return(NULL)
  }

  x_treat <- x_mat_full[idx_treat, , drop = FALSE]
  x_control <- x_mat_full[idx_control, , drop = FALSE]
  arm_weights <- extremeqr_split_weights(weights, d)

  beta_tau_treat <- extremalqr_fit_beta(y[idx_treat], x_treat, tau, weights = arm_weights$treat)
  beta_tau_control <- extremalqr_fit_beta(y[idx_control], x_control, tau, weights = arm_weights$control)
  beta_mtau_treat <- extremalqr_fit_beta(y[idx_treat], x_treat, mtau, weights = arm_weights$treat)
  beta_mtau_control <- extremalqr_fit_beta(y[idx_control], x_control, mtau, weights = arm_weights$control)
  if (is.null(beta_tau_treat) || is.null(beta_tau_control) ||
      is.null(beta_mtau_treat) || is.null(beta_mtau_control)) {
    return(NULL)
  }

  mean_x_treat <- extremalqr_weighted_mean(x_treat, weights = arm_weights$treat)
  mean_x_control <- extremalqr_weighted_mean(x_control, weights = arm_weights$control)
  spacing_hat <- drop(mean_x_treat %*% (beta_mtau_treat - beta_tau_treat)) -
    drop(mean_x_control %*% (beta_mtau_control - beta_tau_control))
  a_hat <- sqrt(tau * n) / spacing_hat
  if (!is.finite(a_hat) || abs(a_hat) < 1e-8) {
    return(NULL)
  }

  tau_b <- min((tau * n) / b, 0.2)
  mtau_b <- tau_b + (d_coef + spacing_p) / b
  if (!is.finite(tau_b) || tau_b <= 0 || tau_b >= 1 || mtau_b >= 1) {
    return(NULL)
  }

  beta_tau_b_treat_full <- if (abs(tau_b - tau) < 1e-10) {
    beta_tau_treat
  } else {
    extremalqr_fit_beta(y[idx_treat], x_treat, tau_b, weights = arm_weights$treat)
  }
  beta_tau_b_control_full <- if (abs(tau_b - tau) < 1e-10) {
    beta_tau_control
  } else {
    extremalqr_fit_beta(y[idx_control], x_control, tau_b, weights = arm_weights$control)
  }
  if (is.null(beta_tau_b_treat_full) || is.null(beta_tau_b_control_full)) {
    return(NULL)
  }

  theta_hat <- as.numeric(phi_mat %*% beta_tau_treat) - as.numeric(phi_mat %*% beta_tau_control)
  theta_tau_b_full <- as.numeric(phi_mat %*% beta_tau_b_treat_full) -
    as.numeric(phi_mat %*% beta_tau_b_control_full)

  subsample_stats <- matrix(NA_real_, nrow = nrow(phi_mat), ncol = n_subsamples)

  for (s in seq_len(n_subsamples)) {
    idx <- sample.int(n, size = b, replace = FALSE)
    y_sub <- y[idx]
    d_sub <- d[idx]
    x_sub_df <- x_df[idx, , drop = FALSE]

    if (length(unique(d_sub)) < 2L) {
      next
    }

    weights_sub <- NULL
    if (reestimate_ps) {
      weights_sub <- extremalqr_compute_ps_weights(
        d = d_sub,
        x = x_sub_df,
        trim = ps_trim,
        cap = ps_cap
      )
    } else if (!is.null(weights)) {
      weights_sub <- weights[idx]
    }
    arm_weights_sub <- extremeqr_split_weights(weights_sub, d_sub)

    idx_treat_sub <- d_sub == 1
    idx_control_sub <- d_sub == 0
    if (!any(idx_treat_sub) || !any(idx_control_sub)) {
      next
    }

    x_treat_sub <- tryCatch(
      extremeqr_model_matrix(prep$terms_x, x_sub_df[idx_treat_sub, , drop = FALSE]),
      error = function(e) NULL
    )
    x_control_sub <- tryCatch(
      extremeqr_model_matrix(prep$terms_x, x_sub_df[idx_control_sub, , drop = FALSE]),
      error = function(e) NULL
    )
    if (is.null(x_treat_sub) || is.null(x_control_sub)) {
      next
    }

    beta_tau_treat_sub <- extremalqr_fit_beta(
      y_sub[idx_treat_sub], x_treat_sub, tau_b, weights = arm_weights_sub$treat
    )
    beta_tau_control_sub <- extremalqr_fit_beta(
      y_sub[idx_control_sub], x_control_sub, tau_b, weights = arm_weights_sub$control
    )
    beta_mtau_treat_sub <- extremalqr_fit_beta(
      y_sub[idx_treat_sub], x_treat_sub, mtau_b, weights = arm_weights_sub$treat
    )
    beta_mtau_control_sub <- extremalqr_fit_beta(
      y_sub[idx_control_sub], x_control_sub, mtau_b, weights = arm_weights_sub$control
    )
    if (is.null(beta_tau_treat_sub) || is.null(beta_tau_control_sub) ||
        is.null(beta_mtau_treat_sub) || is.null(beta_mtau_control_sub)) {
      next
    }

    mean_x_treat_sub <- extremalqr_weighted_mean(x_treat_sub, weights = arm_weights_sub$treat)
    mean_x_control_sub <- extremalqr_weighted_mean(x_control_sub, weights = arm_weights_sub$control)
    spacing_sub <- drop(mean_x_treat_sub %*% (beta_mtau_treat_sub - beta_tau_treat_sub)) -
      drop(mean_x_control_sub %*% (beta_mtau_control_sub - beta_tau_control_sub))
    a_sub <- sqrt(tau_b * b) / spacing_sub
    if (!is.finite(a_sub) || abs(a_sub) < 1e-8) {
      next
    }

    theta_sub <- as.numeric(phi_mat %*% beta_tau_treat_sub) -
      as.numeric(phi_mat %*% beta_tau_control_sub)
    subsample_stats[, s] <- as.numeric(a_sub) * (theta_sub - theta_tau_b_full)
  }

  valid_cols <- colSums(is.finite(subsample_stats)) == nrow(subsample_stats)
  subsample_stats <- subsample_stats[, valid_cols, drop = FALSE]
  if (ncol(subsample_stats) < max(25, ceiling(n_subsamples / 5))) {
    return(NULL)
  }

  c_med <- apply(subsample_stats, 1, stats::quantile, probs = 0.5, na.rm = TRUE, names = FALSE)
  c_lo <- apply(subsample_stats, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
  c_hi <- apply(subsample_stats, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  sd_hat <- apply(subsample_stats, 1, stats::sd, na.rm = TRUE) / abs(a_hat)

  est <- theta_hat - c_med / a_hat
  lo_raw <- theta_hat - c_hi / a_hat
  hi_raw <- theta_hat - c_lo / a_hat

  data.frame(
    estimate = est,
    sd = sd_hat,
    lower = pmin(lo_raw, hi_raw),
    upper = pmax(lo_raw, hi_raw)
  )
}

extremeqr_arm_cqte_core <- function(y, d, x, new_x, tau, alpha = 0.1,
                                    central_boot = 100, extremal_threshold = 20,
                                    spacing_p = 5, b = NULL,
                                    n_subsamples = 200, seed = NULL,
                                    weights = NULL, reestimate_ps = FALSE,
                                    ps_trim = 1e-4, ps_cap = NULL) {
  prep <- extremeqr_prepare_common_model(x)
  x_mat_full <- extremeqr_model_matrix(prep$terms_x, prep$x_df)
  d_coef <- ncol(x_mat_full)

  idx_treat <- d == 1
  idx_control <- d == 0
  use_extremal <- extremalqr_should_use_extremal(
    tau = tau, n = sum(idx_treat), d = d_coef, threshold = extremal_threshold
  ) && extremalqr_should_use_extremal(
    tau = tau, n = sum(idx_control), d = d_coef, threshold = extremal_threshold
  )
  if (!use_extremal) {
    return(NULL)
  }

  if (tau > 0.5) {
    fit <- extremeqr_arm_sn_interval(
      y = -y,
      d = d,
      x = x,
      new_x = new_x,
      tau = 1 - tau,
      alpha = alpha,
      spacing_p = spacing_p,
      b = b,
      n_subsamples = n_subsamples,
      seed = seed,
      weights = weights,
      reestimate_ps = reestimate_ps,
      ps_trim = ps_trim,
      ps_cap = ps_cap
    )
    if (is.null(fit)) {
      return(NULL)
    }
    return(data.frame(
      estimate = -fit$estimate,
      sd = fit$sd,
      lower = -fit$upper,
      upper = -fit$lower
    ))
  }

  extremeqr_arm_sn_interval(
    y = y,
    d = d,
    x = x,
    new_x = new_x,
    tau = tau,
    alpha = alpha,
    spacing_p = spacing_p,
    b = b,
    n_subsamples = n_subsamples,
    seed = seed,
    weights = weights,
    reestimate_ps = reestimate_ps,
    ps_trim = ps_trim,
    ps_cap = ps_cap
  )
}
