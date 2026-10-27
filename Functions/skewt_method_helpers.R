# Log-skew-t regression helpers for the Skew-t competitor (qte_skewt_multi() in Functions/other_method.R).
skewt_require_package <- function() {
  if (!requireNamespace("sn", quietly = TRUE)) {
    stop("Package 'sn' is required for the SkewT method.")
  }
}

skewt_validate_positive_response <- function(y, method_label = "SkewT") {
  if (any(!is.finite(y))) {
    stop(method_label, " requires finite outcome values.")
  }
  if (any(y <= 0)) {
    stop(method_label, " requires strictly positive outcomes because the paper specifies a log-skew-t model.")
  }
}

skewt_make_formula <- function(response_name, x_names) {
  rhs <- if (length(x_names)) paste(x_names, collapse = " + ") else "1"
  stats::as.formula(paste(response_name, "~", rhs))
}

skewt_fit_arm <- function(y, x_df, weights = NULL) {
  skewt_require_package()

  x_df <- as.data.frame(x_df, check.names = FALSE)
  if (is.null(names(x_df))) {
    names(x_df) <- paste0("X", seq_len(ncol(x_df)))
  }

  fit_df <- data.frame(y = y, x_df, check.names = FALSE)
  fit_formula <- skewt_make_formula("y", names(x_df))
  fit_weights <- NULL
  if (!is.null(weights)) {
    fit_weights <- as.numeric(weights)
    if (length(fit_weights) != nrow(fit_df) || any(!is.finite(fit_weights))) {
      return(NULL)
    }
    fit_df$.weights <- fit_weights
  }

  # a rare bootstrap sample makes selm run for a very long time: give up after 20 seconds
  setTimeLimit(elapsed = 20, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  suppressWarnings(tryCatch(
    if (is.null(fit_weights)) {
      sn::selm(fit_formula, family = "ST", data = fit_df)
    } else {
      sn::selm(fit_formula, family = "ST", data = fit_df, weights = .weights)
    },
    error = function(e) NULL
  ))
}

skewt_error_dp <- function(fit) {
  dp <- tryCatch(stats::coef(fit, param.type = "DP"), error = function(e) NULL)
  if (is.null(dp)) {
    return(NULL)
  }

  nms <- names(dp)
  omega <- dp[grepl("^omega", nms)]
  alpha <- dp[grepl("^alpha", nms)]
  nu <- dp[grepl("^nu", nms)]

  if (!length(omega) || !length(alpha) || !length(nu)) {
    return(NULL)
  }

  c(0, unname(omega[[1]]), unname(alpha[[1]]), unname(nu[[1]]))
}

skewt_predict_quantiles <- function(fit, new_x, tau_vector) {
  skewt_require_package()
  setTimeLimit(elapsed = 20, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)

  tau_vector <- as.numeric(tau_vector)
  new_x <- as.data.frame(new_x, check.names = FALSE)
  if (is.null(names(new_x))) {
    names(new_x) <- paste0("X", seq_len(ncol(new_x)))
  }

  fit_terms <- tryCatch(
    stats::delete.response(attr(fit@input$model, "terms")),
    error = function(e) NULL
  )
  dp <- tryCatch(
    stats::coef(fit, param.type = "DP"),
    error = function(e) NULL
  )
  linear_pred <- NULL
  if (!is.null(fit_terms) && !is.null(dp)) {
    x_mat <- tryCatch(
      stats::model.matrix(fit_terms, data = new_x),
      error = function(e) NULL
    )
    beta <- dp[!grepl("^(omega|alpha|nu)", names(dp))]
    beta_names <- names(beta)
    beta_names[beta_names == "(Intercept.DP)"] <- "(Intercept)"
    beta_names <- sub("\\.DP$", "", beta_names)
    if (!is.null(x_mat) && length(beta) && all(colnames(x_mat) %in% beta_names)) {
      beta_vec <- as.numeric(beta[match(colnames(x_mat), beta_names)])
      linear_pred <- as.numeric(x_mat %*% beta_vec)
    }
  }
  error_dp <- skewt_error_dp(fit)

  if (is.null(linear_pred) || is.null(error_dp)) {
    return(NULL)
  }

  # Error quantiles by Monte Carlo: sn::qst can hang (in compiled code, so it cannot be interrupted)
  # at extreme tau for large shape and degrees of freedom; 2e6 draws leave about 180 above tau = 0.99991.
  error_dp[[4]] <- min(error_dp[[4]], 1e4)
  error_quantiles <- tryCatch(
    as.numeric(stats::quantile(sn::rst(2e6, dp = error_dp), tau_vector, names = FALSE)),
    error = function(e) NULL
  )
  if (is.null(error_quantiles)) {
    return(NULL)
  }

  preds <- outer(linear_pred, error_quantiles, "+")
  dimnames(preds) <- list(NULL, NULL)
  preds
}
