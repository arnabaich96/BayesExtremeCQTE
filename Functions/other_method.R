# Competing methods of the simulation study (Simulation/Competing.R):
#   qte_ipw()         IPW linear quantile regression with treatment-by-covariate interactions
#   qte_extremeqr()   weighted extreme quantile regression (Functions/extremal_qr_inference.R)
#   qte_skewt_multi() arm-wise log-skew-t regression, bootstrap intervals (Functions/skewt_method_helpers.R)
# Each returns, per covariate profile, the estimate, interval, bias, coverage and interval length.
source(here::here("Functions/extremal_qr_inference.R"))
source(here::here("Functions/skewt_method_helpers.R"))

summarize_boot_percentile <- function(boot_draws, true_qte, alpha = 0.05, probs = NULL,
                                      point_estimate = NULL) {
  if (is.null(probs)) {
    probs <- rep(0.5, nrow(boot_draws))
  }
  if (length(probs) == 1L) {
    probs <- rep(probs, nrow(boot_draws))
  }

  if (is.null(point_estimate)) {
    estimate <- vapply(seq_len(nrow(boot_draws)), function(i) {
      stats::quantile(boot_draws[i, ], probs = probs[[i]], na.rm = TRUE, names = FALSE)
    }, numeric(1))
  } else {
    estimate <- as.numeric(point_estimate)
    if (length(estimate) != nrow(boot_draws)) {
      stop("point_estimate must have length equal to nrow(boot_draws).")
    }
  }
  lower <- apply(boot_draws, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
  upper <- apply(boot_draws, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  sd_hat <- apply(boot_draws, 1, stats::sd, na.rm = TRUE)

  round(data.frame(
    estimate = estimate,
    lower = lower,
    upper = upper,
    bias = estimate - true_qte,
    coverage = as.numeric(lower <= true_qte & true_qte <= upper),
    length_ci = upper - lower,
    sd = sd_hat
  ), 3)
}

qte_extremeqr <- function(Y, D, X, new_X, tau, B = 200, true_qte, alpha) {
  weights <- extremalqr_compute_ps_weights(d = D, x = X)
  fit <- extremeqr_arm_cqte_core(
    y = Y,
    d = D,
    x = X,
    new_x = new_X,
    tau = tau,
    alpha = alpha,
    central_boot = max(50, floor(B / 2)),
    n_subsamples = B,
    weights = weights,
    reestimate_ps = TRUE
  )

  if (is.null(fit)) {
    fit <- extremeqr_arm_boot_interval(
      y = Y,
      d = D,
      x = X,
      new_x = new_X,
      tau = tau,
      alpha = alpha,
      n_boot = max(50, floor(B / 2)),
      weights = weights,
      reestimate_ps = TRUE
    )
    if (is.null(fit)) {
      return(qte_ipw(
        Y = Y,
        D = D,
        X = X,
        new_X = new_X,
        tau = tau,
        true_qte = true_qte,
        alpha = alpha
      ))
    }
  }

  round(data.frame(
    estimate = fit$estimate,
    lower = fit$lower,
    upper = fit$upper,
    bias = fit$estimate - true_qte,
    coverage = as.numeric(fit$lower <= true_qte & true_qte <= fit$upper),
    length_ci = fit$upper - fit$lower
  ), 3)
}

# Same estimator as qte_skewt(), for several tau at once: the arm-wise fits (original sample and B
# bootstrap samples) do not depend on tau, so they are computed once. true_qte: G x length(taus).
# Returns a list over taus of the data frames returned by qte_skewt().
qte_skewt_multi <- function(Y, D, X, new_X, taus, B = 100, true_qte, alpha) {
  X <- as.data.frame(X, check.names = FALSE)
  new_X <- as.data.frame(new_X, check.names = FALSE)
  colnames(new_X) <- colnames(X)
  skewt_validate_positive_response(Y, method_label = "SkewT")
  G <- nrow(new_X); n_tau <- length(taus)
  na_list <- function() {
    na <- rep(NA_real_, G)
    rep(list(data.frame(estimate = na, lower = na, upper = na, bias = na, coverage = na, length_ci = na)), n_tau)
  }
  cqte_all <- function(y, d, x) tryCatch({  # G x n_tau matrix of CQTE estimates, or NULL
    f1 <- skewt_fit_arm(log(y[d == 1]), x[d == 1, , drop = FALSE])
    f0 <- skewt_fit_arm(log(y[d == 0]), x[d == 0, , drop = FALSE])
    if (is.null(f1) || is.null(f0)) return(NULL)
    q1 <- skewt_predict_quantiles(f1, new_X, taus)
    q0 <- skewt_predict_quantiles(f0, new_X, taus)
    if (is.null(q1) || is.null(q0)) return(NULL)
    exp(q1) - exp(q0)
  }, error = function(e) NULL)              # includes the 20-second time limits of the helpers
  point <- cqte_all(Y, D, X)
  if (is.null(point)) return(na_list())

  boots <- array(NA_real_, dim = c(G, n_tau, B))
  for (b in seq_len(B)) {
    idx <- sample.int(length(Y), replace = TRUE)
    if (sum(D[idx] == 1) <= 1 || sum(D[idx] == 0) <= 1) next
    est <- cqte_all(Y[idx], D[idx], X[idx, , drop = FALSE])
    if (!is.null(est)) boots[, , b] <- est
  }

  lapply(seq_len(n_tau), function(t) {
    bt <- boots[, t, , drop = TRUE]
    bt <- bt[, colSums(is.finite(bt)) == G, drop = FALSE]
    if (ncol(bt) < max(25, ceiling(B / 5))) return(na_list()[[1]])
    s <- summarize_boot_percentile(bt, true_qte = true_qte[, t], alpha = alpha, probs = taus[t],
                                   point_estimate = point[, t])
    s[, c("estimate", "lower", "upper", "bias", "coverage", "length_ci")]
  })
}

qte_ipw <- function(Y, D, X, new_X, tau, true_qte, alpha) {
  data <- data.frame(Y = Y, D = D, X)
  ps_model <- glm(D ~ ., data = data, family = binomial)
  ps <- predict(ps_model, type = "response")
  weights <- ifelse(D == 1, 1 / ps, 1 / (1 - ps))
  
  # D interacts with every covariate, so the CQTE at x is the contrast D + sum_k x_k D:X_k
  formula <- as.formula(paste("Y ~ D * (", paste(names(X), collapse = " + "), ")"))
  model <- rq(formula, tau = tau, data = data, weights = weights)
  pred_treated <- predict(model, newdata = cbind(D = 1, new_X))
  pred_control <- predict(model, newdata = cbind(D = 0, new_X))
  qte <- pred_treated - pred_control

  # standard error of the contrast from the bootstrap covariance of the coefficients
  V <- summary(model, se = "boot", covariance = TRUE)$cov
  C <- matrix(0, nrow(new_X), length(coef(model)), dimnames = list(NULL, names(coef(model))))
  C[, "D"] <- 1
  for (v in names(X)) C[, paste0("D:", v)] <- new_X[[v]]
  se <- sqrt(rowSums((C %*% V) * C))
  ci_lower <- qte - qnorm(1 - alpha / 2) * se
  ci_upper <- qte + qnorm(1 - alpha / 2) * se
  bias <- qte - true_qte
  coverage <- as.numeric(ci_lower <= true_qte & true_qte <= ci_upper)
  length_ci <- ci_upper - ci_lower
  
  round(data.frame(
    estimate = qte,
    lower = ci_lower,
    upper = ci_upper,
    bias = bias,
    coverage = coverage,
    length_ci = length_ci
  ), 3)
}






